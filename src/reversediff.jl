#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#
# Direct ReverseDiff dispatch for FMIBase.eval!
# Bypasses ChainRules rrule for better performance and fewer allocations.

# Cache struct for direct ReverseDiff dispatch.
# Captures primal state at forward pass time + pre-allocated gradient buffers.
mutable struct FMIEvalRDCache
    c::Any              # FMUInstance
    cRef::UInt64
    # Primal values at forward-pass time (copies for correct pullback)
    x_val::Vector{Float64}
    u_val::Vector{Float64}
    p_val::Vector{Float64}
    t_val::Float64
    x_d::Vector{Float64}
    # Reference arrays
    u_refs::Vector{UInt32}
    p_refs::Vector{UInt32}
    dx_refs::Vector{UInt32}
    y_refs::Vector{UInt32}
    ec_idcs::Vector{UInt32}
    # Lengths for output slicing
    dx_len::Int
    y_len::Int
    ec_len::Int
    # Flags for which inputs are tracked
    has_x::Bool
    has_u::Bool
    has_p::Bool
    has_t::Bool
    # Pre-allocated gradient accumulators (zero alloc in backward pass)
    x̄::Vector{Float64}
    ū::Vector{Float64}
    p̄::Vector{Float64}
    t̄::Vector{Float64}
    # Pre-allocated back_output buffers (avoids zeros() per backward pass)
    back_dx::Vector{Float64}
    back_y::Vector{Float64}
    back_ec::Vector{Float64}
    back_xd::Vector{Float64}
end

"""
    _eval_reversediff!(c, cRef, dx, dx_refs, y, y_refs, x_val, u_val, u_refs, p_val, p_refs,
                        ec_val, ec_idcs, t_val, x_d, tracked_inputs, tp,
                        has_x, has_u, has_p, has_t)

Shared kernel for all ReverseDiff dispatch methods.
"""
function _eval_reversediff!(
    c, cRef::UInt64,
    dx::AbstractVector{Float64}, dx_refs::Vector{UInt32},
    y::AbstractVector{Float64}, y_refs::Vector{UInt32},
    x_val::Vector{Float64}, u_val::Vector{Float64}, u_refs::Vector{UInt32},
    p_val::Vector{Float64}, p_refs::Vector{UInt32},
    ec_val::Vector{Float64}, ec_idcs::Vector{UInt32},
    t_val::Float64, x_d::Vector{Float64},
    tracked_inputs::Tuple, tp,
    has_x::Bool, has_u::Bool, has_p::Bool, has_t::Bool,
)
    # --- Primal evaluation ---
    Ω = FMIBase.eval!(
        cRef, dx, dx_refs, y, y_refs, x_val, u_val, u_refs,
        p_val, p_refs, ec_val, ec_idcs, t_val, x_d,
    )

    # Track output on the tape
    output = ReverseDiff.track(Ω, tp)

    dx_len = length(dx)
    y_len = length(y)
    ec_len = length(ec_val)

    # Build cache with copies of primal state and pre-allocated gradient buffers
    cache = FMIEvalRDCache(
        c, cRef,
        copy(x_val), copy(u_val), copy(p_val), t_val, copy(x_d),
        copy(u_refs), copy(p_refs), copy(dx_refs), copy(y_refs), copy(ec_idcs),
        dx_len, y_len, ec_len,
        has_x, has_u, has_p, has_t,
        zeros(Float64, length(x_val)),  # x̄
        zeros(Float64, length(u_val)),  # ū
        zeros(Float64, length(p_val)),  # p̄
        zeros(Float64, 1),              # t̄
        zeros(Float64, dx_len),         # back_dx
        zeros(Float64, y_len),          # back_y
        zeros(Float64, ec_len),         # back_ec
        zeros(Float64, length(x_d)),    # back_xd
    )

    ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, FMIBase.eval!, tracked_inputs, output, cache)
    return output
end

# Backward pass: dispatch on our FMIEvalRDCache cache type
@noinline function ReverseDiff.special_reverse_exec!(
    instruction::ReverseDiff.SpecialInstruction{typeof(FMIBase.eval!),<:Tuple,<:Any,FMIEvalRDCache},
)
    pb = instruction.cache
    c = pb.c
    r̄ = ReverseDiff.deriv(instruction.output)

    # Slice output adjoint into dx, y, ec portions
    d̄x = @view(r̄[1:pb.dx_len])
    ȳ = @view(r̄[(pb.dx_len+1):(pb.dx_len+pb.y_len)])
    ēc = @view(r̄[(pb.dx_len+pb.y_len+1):end])

    # Reset gradient accumulators
    fill!(pb.x̄, 0.0)
    fill!(pb.ū, 0.0)
    fill!(pb.p̄, 0.0)
    fill!(pb.t̄, 0.0)

    # Restore FMU state to forward-pass state
    FMIBase.eval_set!(c, pb.x_val, pb.u_val, pb.u_refs, pb.p_val, pb.p_refs, pb.t_val, pb.x_d)

    x_refs = c.fmu.modelDescription.stateValueReferences
    dx_refs = pb.dx_refs
    y_refs = pb.y_refs
    ec_idcs = pb.ec_idcs
    u_refs = pb.u_refs
    p_refs = pb.p_refs

    states = length(pb.x_val) > 0
    inputs = length(u_refs) > 0
    parameters = length(p_refs) > 0
    times = FMIBase.isSetReal(c.fmu, pb.t_val)

    derivatives = pb.dx_len > 0 && !isZeroTangent(d̄x) && any(!iszero, d̄x)
    outputs = pb.y_len > 0 && !isZeroTangent(ȳ) && any(!iszero, ȳ)
    eventIndicators = pb.ec_len > 0 && length(ec_idcs) > 0

    # Resolve dx_refs
    if pb.dx_len > 0 && length(dx_refs) == 0
        dx_refs = c.fmu.modelDescription.derivativeValueReferences
    end

    # VJP computations
    if derivatives
        if pb.has_x && states
            vjp!(c, :∂ẋ_∂x, dx_refs, x_refs, pb.x_val, d̄x; accu = pb.x̄)
            c.solution.evals_∂ẋ_∂x += 1
        end
        if pb.has_u && inputs
            vjp!(c, :∂ẋ_∂u, dx_refs, u_refs, pb.u_val, d̄x; accu = pb.ū)
            c.solution.evals_∂ẋ_∂u += 1
        end
        if pb.has_p && parameters
            vjp!(c, :∂ẋ_∂p, dx_refs, p_refs, pb.p_val, d̄x; accu = pb.p̄)
            c.solution.evals_∂ẋ_∂p += 1
        end
        if pb.has_t && times && c.fmu.executionConfig.eval_t_gradients
            vgp!(c, :∂ẋ_∂t, dx_refs, :time, pb.t_val, d̄x; accu = pb.t̄)
            c.solution.evals_∂ẋ_∂t += 1
        end
    end

    if outputs
        if pb.has_x && states
            vjp!(c, :∂y_∂x, y_refs, x_refs, pb.x_val, ȳ; accu = pb.x̄)
            c.solution.evals_∂y_∂x += 1
        end
        if pb.has_u && inputs
            vjp!(c, :∂y_∂u, y_refs, u_refs, pb.u_val, ȳ; accu = pb.ū)
            c.solution.evals_∂y_∂u += 1
        end
        if pb.has_p && parameters
            vjp!(c, :∂y_∂p, y_refs, p_refs, pb.p_val, ȳ; accu = pb.p̄)
            c.solution.evals_∂y_∂p += 1
        end
        if pb.has_t && times && c.fmu.executionConfig.eval_t_gradients
            vgp!(c, :∂y_∂t, y_refs, :time, pb.t_val, ȳ; accu = pb.t̄)
            c.solution.evals_∂y_∂t += 1
        end
    end

    if eventIndicators
        if pb.has_x && states
            vjp!(c, :∂e_∂x, (:indicators, ec_idcs), x_refs, pb.x_val, ēc; accu = pb.x̄)
            c.solution.evals_∂e_∂x += 1
        end
        if pb.has_u && inputs
            vjp!(c, :∂e_∂u, (:indicators, ec_idcs), u_refs, pb.u_val, ēc; accu = pb.ū)
            c.solution.evals_∂e_∂u += 1
        end
        if pb.has_p && parameters
            vjp!(c, :∂e_∂p, (:indicators, ec_idcs), p_refs, pb.p_val, ēc; accu = pb.p̄)
            c.solution.evals_∂e_∂p += 1
        end
        if pb.has_t && times && c.fmu.executionConfig.eval_t_gradients
            vgp!(c, :∂e_∂t, (:indicators, ec_idcs), :time, pb.t_val, ēc; accu = pb.t̄)
            c.solution.evals_∂e_∂t += 1
        end
    end

    # Distribute gradients back to tracked inputs via the rrule return convention.
    # Reuse pre-allocated buffers for zero-gradient outputs (avoids zeros() per backward pass)
    fill!(pb.back_dx, 0.0)
    fill!(pb.back_y, 0.0)
    fill!(pb.back_ec, 0.0)
    fill!(pb.back_xd, 0.0)

    # Input derivs matching: (c̄Ref, d̄x, d̄x_refs, ȳ, ȳ_refs, x̄, ū, ū_refs, p̄, p̄_refs, ēc, ēc_idcs, t̄, x̄_d)
    input_derivs = (
        [],                     # c̄Ref
        pb.back_dx,             # d̄x (output, not differentiated w.r.t.)
        [],                     # d̄x_refs
        pb.back_y,              # ȳ (output, not differentiated w.r.t.)
        [],                     # ȳ_refs
        pb.x̄,                   # x̄
        pb.ū,                   # ū
        [],                     # ū_refs
        pb.p̄,                   # p̄
        [],                     # p̄_refs
        pb.back_ec,             # ēc (output)
        [],                     # ēc_idcs
        pb.t̄[1],                # t̄
        pb.back_xd,             # x̄_d
    )

    ReverseDiff._add_to_deriv!.(instruction.input, input_derivs)
    ReverseDiff.unseed!(instruction.output)
    return nothing
end

# Forward pass for compiled tape re-evaluation
@noinline function ReverseDiff.special_forward_exec!(
    instruction::ReverseDiff.SpecialInstruction{typeof(FMIBase.eval!),<:Tuple,<:Any,FMIEvalRDCache},
)
    output, input = instruction.output, instruction.input
    ReverseDiff.pull_value!.(input)

    pb = instruction.cache

    # Re-extract current values from tracked inputs
    input_vals = map(ReverseDiff.value, input)

    # Re-run the primal evaluation
    Ω = FMIBase.eval!(
        pb.cRef,
        zeros(Float64, pb.dx_len), pb.dx_refs,
        zeros(Float64, pb.y_len), pb.y_refs,
        Float64.(input_vals[6]),  # x
        Float64.(input_vals[7]),  # u
        pb.u_refs,
        Float64.(input_vals[9]),  # p
        pb.p_refs,
        zeros(Float64, pb.ec_len), pb.ec_idcs,
        Float64(input_vals[13]),  # t
        pb.x_d,
    )
    ReverseDiff.value!(output, Ω)
    return nothing
end

# ============================================================================
# Entry-point methods for each argument type combination.
# ============================================================================

# Macro for ReverseDiff dispatch methods.
# Takes 7 args: (dx_type, y_type, x_type, u_type, p_type, ec_type, t_type)
# matching the original @grad_from_chainrules signatures.
macro _rd_eval(dx_type, y_type, x_type, u_type, p_type, ec_type, t_type)
    dx_T = dx_type == :Tracked ? :(AbstractVector{<:ReverseDiff.TrackedReal}) : :(AbstractVector{<:Real})
    y_T = y_type == :Tracked ? :(AbstractVector{<:ReverseDiff.TrackedReal}) : :(AbstractVector{<:Real})
    x_T = x_type == :Tracked ? :(AbstractVector{<:ReverseDiff.TrackedReal}) : :(AbstractVector{<:Real})
    u_T = u_type == :Tracked ? :(AbstractVector{<:ReverseDiff.TrackedReal}) : :(AbstractVector{<:Real})
    p_T = p_type == :Tracked ? :(AbstractVector{<:ReverseDiff.TrackedReal}) : :(AbstractVector{<:Real})
    ec_T = ec_type == :Tracked ? :(AbstractVector{<:ReverseDiff.TrackedReal}) : :(AbstractVector{<:Real})
    t_T = t_type == :Tracked ? :(ReverseDiff.TrackedReal) : :Real

    # Note: y_refs uses UInt32 for ReverseDiff (matching original @grad_from_chainrules)
    y_refs_T = :(AbstractVector{<:UInt32})
    u_refs_T = :(AbstractVector{<:UInt32})
    p_refs_T = :(AbstractVector{<:UInt32})

    has_x = x_type == :Tracked
    has_u = u_type == :Tracked
    has_p = p_type == :Tracked
    has_t = t_type == :Tracked

    # Value extraction expressions (avoid allocating copies when possible)
    x_val = x_type == :Tracked ? :(_rd_extract_values(x)) : :(_as_f64_vec(x))
    u_val = u_type == :Tracked ? :(_rd_extract_values(u)) : :(_as_f64_vec(u))
    p_val = p_type == :Tracked ? :(_rd_extract_values(p)) : :(_as_f64_vec(p))
    ec_val = ec_type == :Tracked ? :(_rd_extract_values(ec)) : :(_as_f64_vec(ec))
    t_val_expr = t_type == :Tracked ? :(Float64(ReverseDiff.value(t))) : :(Float64(t))

    # Collect all tracked args for tape access
    tracked_vars = []
    x_type == :Tracked && push!(tracked_vars, :x)
    u_type == :Tracked && push!(tracked_vars, :u)
    p_type == :Tracked && push!(tracked_vars, :p)
    ec_type == :Tracked && push!(tracked_vars, :ec)
    t_type == :Tracked && push!(tracked_vars, :t)
    dx_type == :Tracked && push!(tracked_vars, :dx)
    y_type == :Tracked && push!(tracked_vars, :y)

    # Get tape from first tracked variable
    tape_expr = :(ReverseDiff.tape($(tracked_vars[1])...))

    # For dx/y: if tracked, allocate Float64 buffers for primal
    dx_buf_expr = dx_type == :Tracked ? :(zeros(Float64, length(dx))) : :dx
    y_buf_expr = y_type == :Tracked ? :(zeros(Float64, length(y))) : :y

    quote
        function FMIBase.eval!(
            cRef::UInt64,
            dx::$dx_T,
            dx_refs::AbstractVector{<:fmiValueReference},
            y::$y_T,
            y_refs::$y_refs_T,
            x::$x_T,
            u::$u_T,
            u_refs::$u_refs_T,
            p::$p_T,
            p_refs::$p_refs_T,
            ec::$ec_T,
            ec_idcs::AbstractVector{<:fmiValueReference},
            t::$t_T,
            x_d::AbstractVector{<:Real},
        )
            c = unsafe_pointer_to_objref(Ptr{Nothing}(cRef))
            tp = $tape_expr

            x_val = $x_val
            u_val = $u_val
            p_val = $p_val
            ec_val = $ec_val
            t_val = $t_val_expr
            x_d_clean = _as_f64_vec(x_d)

            dx_refs_clean = _as_u32_vec(dx_refs)
            y_refs_clean = _as_u32_vec(y_refs)
            u_refs_clean = _as_u32_vec(u_refs)
            p_refs_clean = _as_u32_vec(p_refs)
            ec_idcs_clean = _as_u32_vec(ec_idcs)

            dx_buf = $dx_buf_expr
            y_buf = $y_buf_expr

            # Build the tracked_inputs tuple matching the full eval! signature
            # (cRef, dx, dx_refs, y, y_refs, x, u, u_refs, p, p_refs, ec, ec_idcs, t, x_d)
            tracked_inputs = (cRef, dx, dx_refs, y, y_refs, x, u, u_refs, p, p_refs, ec, ec_idcs, t, x_d)

            return _eval_reversediff!(
                c, cRef,
                dx_buf, dx_refs_clean, y_buf, y_refs_clean,
                x_val, u_val, u_refs_clean, p_val, p_refs_clean,
                ec_val, ec_idcs_clean, t_val, x_d_clean,
                tracked_inputs, tp,
                $has_x, $has_u, $has_p, $has_t,
            )
        end
    end |> esc
end

# All 17 combinations matching the original @grad_from_chainrules registrations exactly.
# Format: @_rd_eval dx y x u p ec t

# 1. dx, y, x, u, p, ec, t (all tracked)
@_rd_eval Tracked Tracked Tracked Tracked Tracked Tracked Tracked

# 2. dx, y, x, u, t
@_rd_eval Tracked Tracked Tracked Tracked Real Real Tracked

# 3. x, u
@_rd_eval Real Real Tracked Tracked Real Real Real

# 4. x, u, t
@_rd_eval Real Real Tracked Tracked Real Real Tracked

# 5. x, p
@_rd_eval Real Real Tracked Real Tracked Real Real

# 6. t only
@_rd_eval Real Real Real Real Real Real Tracked

# 7. x only
@_rd_eval Real Real Tracked Real Real Real Real

# 8. u only
@_rd_eval Real Real Real Tracked Real Real Real

# 9. p only
@_rd_eval Real Real Real Real Tracked Real Real

# 10. ec only
@_rd_eval Real Real Real Real Real Tracked Real

# 11. x, t
@_rd_eval Real Real Tracked Real Real Real Tracked

# 12. x, ec, t
@_rd_eval Real Real Tracked Real Real Tracked Tracked

# 13. ec, t
@_rd_eval Real Real Real Real Real Tracked Tracked

# 14. x, ec
@_rd_eval Real Real Tracked Real Real Tracked Real

# 15. x, p, t
@_rd_eval Real Real Tracked Real Tracked Real Tracked

# 16. x, p, ec, t
@_rd_eval Real Real Tracked Real Tracked Tracked Tracked

# 17. x, p, ec
@_rd_eval Real Real Tracked Real Tracked Tracked Real
