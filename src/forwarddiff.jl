#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#
# Direct ForwardDiff dispatch for FMIBase.eval!
# Bypasses ChainRules frule for better performance and fewer allocations.

# Pre-allocated cache for ForwardDiff partial accumulation per FMU instance.
mutable struct FMIFDPartialCache
    ∂dx::Vector{Float64}
    ∂y::Vector{Float64}
    ∂e::Vector{Float64}
    x_val::Vector{Float64}
    u_val::Vector{Float64}
    p_val::Vector{Float64}
    ec_val::Vector{Float64}
    # Cached primal copies (avoid copy() per call)
    dx_primal::Vector{Float64}
    y_primal::Vector{Float64}
    ec_primal::Vector{Float64}
    # Cached partial matrices (avoid Matrix alloc per call)
    ∂dx_parts::Matrix{Float64}
    ∂y_parts::Matrix{Float64}
    ∂e_parts::Matrix{Float64}
    last_N::Int  # track last N to know when to reallocate matrices
end

const _fd_cache = WeakKeyDict{Any,FMIFDPartialCache}()
const _fd_cache_lock = ReentrantLock()

function _get_fd_cache(c)
    lock(_fd_cache_lock) do
        get!(_fd_cache, c) do
            FMIFDPartialCache(
                Float64[], Float64[], Float64[],
                Float64[], Float64[], Float64[], Float64[],
                Float64[], Float64[], Float64[],
                Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0),
                0,
            )
        end
    end
end

# Fill pre-allocated seed vector from i-th partial of a Dual vector (zero allocation)
@inline function _fill_seed!(seed::AbstractVector, duals::AbstractVector{<:ForwardDiff.Dual}, i::Int)
    @inbounds for j in eachindex(duals)
        seed[j] = ForwardDiff.partials(duals[j], i)
    end
end

# Resize a buffer in-place, only allocating when size changes
@inline function _ensure_size!(buf::Vector{Float64}, n::Int)
    if length(buf) != n
        resize!(buf, n)
    end
    buf
end

# Extract primal values into a pre-allocated buffer (avoids ForwardDiff.value.(x) allocation)
@inline function _extract_values!(buf::Vector{Float64}, duals::AbstractVector{<:ForwardDiff.Dual})
    _ensure_size!(buf, length(duals))
    @inbounds for j in eachindex(duals)
        buf[j] = ForwardDiff.value(duals[j])
    end
    buf
end

"""
    _eval_forwarddiff!(c, cRef, dx, dx_refs, y, y_refs, x_val, u_val, u_refs, p_val, p_refs,
                        ec_val, ec_idcs, t_val, x_d,
                        x_duals, u_duals, p_duals, ec_duals, t_dual,
                        Val{TAG}, Val{N}) where {TAG, N}

Shared kernel for all ForwardDiff dispatch methods.
- `x_val`, `u_val`, `p_val`, `ec_val`, `t_val` are Float64 primal values
- `x_duals`, `u_duals`, `p_duals`, `ec_duals` are `nothing` or the original Dual vectors
- `t_dual` is `nothing` or the original Dual scalar
"""
function _eval_forwarddiff!(
    c, cRef::UInt64,
    dx::AbstractVector{Float64}, dx_refs::AbstractVector{<:fmiValueReference},
    y::AbstractVector{Float64}, y_refs::AbstractVector{<:fmiValueReference},
    x_val::AbstractVector{Float64}, u_val::AbstractVector{Float64},
    u_refs::AbstractVector{<:fmiValueReference},
    p_val::AbstractVector{Float64}, p_refs::AbstractVector{<:fmiValueReference},
    ec_val::AbstractVector{Float64}, ec_idcs::AbstractVector{<:fmiValueReference},
    t_val::Float64, x_d::AbstractVector{<:Real},
    x_duals, u_duals, p_duals, ec_duals, t_dual,
    ::Val{TAG}, ::Val{N},
) where {TAG,N}

    # --- Primal evaluation ---
    Ω = FMIBase.eval!(
        cRef, dx, dx_refs, y, y_refs, x_val, u_val, u_refs,
        p_val, p_refs, ec_val, ec_idcs, t_val, x_d,
    )

    ndx = length(Ω.dx)
    ny = length(Ω.y)
    ne = length(Ω.ec)

    # --- Get pre-allocated tangent accumulators ---
    cache = _get_fd_cache(c)
    _ensure_size!(cache.∂dx, ndx)
    _ensure_size!(cache.∂y, ny)
    _ensure_size!(cache.∂e, ne)

    # Capture primal values into cache buffers (eval! wrote into dx/y/ec in-place)
    _ensure_size!(cache.dx_primal, ndx)
    copyto!(cache.dx_primal, Ω.dx)
    _ensure_size!(cache.y_primal, ny)
    copyto!(cache.y_primal, Ω.y)
    _ensure_size!(cache.ec_primal, ne)
    copyto!(cache.ec_primal, Ω.ec)

    # Reuse partial matrices from cache (only reallocate if dimensions changed)
    if cache.last_N != N || size(cache.∂dx_parts, 1) != ndx
        cache.∂dx_parts = Matrix{Float64}(undef, ndx, N)
    end
    if cache.last_N != N || size(cache.∂y_parts, 1) != ny
        cache.∂y_parts = Matrix{Float64}(undef, ny, N)
    end
    if cache.last_N != N || size(cache.∂e_parts, 1) != ne
        cache.∂e_parts = Matrix{Float64}(undef, ne, N)
    end
    cache.last_N = N
    ∂dx_parts = cache.∂dx_parts
    ∂y_parts = cache.∂y_parts
    ∂e_parts = cache.∂e_parts

    # Determine ref arrays
    state_refs = c.fmu.modelDescription.stateValueReferences

    has_x = !isnothing(x_duals) && length(x_duals) > 0
    has_u = !isnothing(u_duals) && length(u_duals) > 0
    has_p = !isnothing(p_duals) && length(p_duals) > 0
    has_ec = !isnothing(ec_duals) && length(ec_duals) > 0
    has_t = !isnothing(t_dual)

    outputs = ny > 0
    derivatives = ndx > 0
    eventIndicators = ne > 0 && length(ec_idcs) > 0
    states = length(x_val) > 0
    inputs = length(u_refs) > 0
    parameters = length(p_refs) > 0
    times = FMIBase.isSetReal(c.fmu, t_val)

    # Resolve dx_refs (same logic as frule)
    if derivatives && length(dx_refs) == 0 &&
       ndx == length(c.fmu.modelDescription.derivativeValueReferences)
        dx_refs = c.fmu.modelDescription.derivativeValueReferences
    end

    # --- Compute JVPs for each partial direction ---
    for i in 1:N
        fill!(cache.∂dx, 0.0)
        fill!(cache.∂y, 0.0)
        fill!(cache.∂e, 0.0)

        # State tangent contributions
        if has_x && states
            jac_x = getfield(c, :∂ẋ_∂x)
            seed = if !isnothing(jac_x)
                _ensure_size!(jac_x.seed, length(x_duals))
                _fill_seed!(jac_x.seed, x_duals, i)
                jac_x.seed
            else
                # First call: jac doesn't exist yet, jvp! will create it.
                # Use a temporary seed (allocated once, then cached via jac.seed on next call).
                [ForwardDiff.partials(x_duals[j], i) for j in eachindex(x_duals)]
            end

            if derivatives
                jvp!(c, :∂ẋ_∂x, dx_refs, state_refs, x_val, seed; accu = cache.∂dx)
                c.solution.evals_∂ẋ_∂x += 1
            end
            if outputs
                jvp!(c, :∂y_∂x, y_refs, state_refs, x_val, seed; accu = cache.∂y)
                c.solution.evals_∂y_∂x += 1
            end
            if eventIndicators
                jvp!(c, :∂e_∂x, (:indicators, ec_idcs), state_refs, x_val, seed; accu = cache.∂e)
                c.solution.evals_∂e_∂x += 1
            end
        end

        # Input tangent contributions
        if has_u && inputs
            jac_u = getfield(c, :∂ẋ_∂u)
            seed = if !isnothing(jac_u)
                _ensure_size!(jac_u.seed, length(u_duals))
                _fill_seed!(jac_u.seed, u_duals, i)
                jac_u.seed
            else
                [ForwardDiff.partials(u_duals[j], i) for j in eachindex(u_duals)]
            end

            if derivatives
                jvp!(c, :∂ẋ_∂u, dx_refs, u_refs, u_val, seed; accu = cache.∂dx)
                c.solution.evals_∂ẋ_∂u += 1
            end
            if outputs
                jvp!(c, :∂y_∂u, y_refs, u_refs, u_val, seed; accu = cache.∂y)
                c.solution.evals_∂y_∂u += 1
            end
            if eventIndicators
                jvp!(c, :∂e_∂u, (:indicators, ec_idcs), u_refs, u_val, seed; accu = cache.∂e)
                c.solution.evals_∂e_∂u += 1
            end
        end

        # Parameter tangent contributions
        if has_p && parameters
            jac_p = getfield(c, :∂ẋ_∂p)
            seed = if !isnothing(jac_p)
                _ensure_size!(jac_p.seed, length(p_duals))
                _fill_seed!(jac_p.seed, p_duals, i)
                jac_p.seed
            else
                [ForwardDiff.partials(p_duals[j], i) for j in eachindex(p_duals)]
            end

            if derivatives
                jvp!(c, :∂ẋ_∂p, dx_refs, p_refs, p_val, seed; accu = cache.∂dx)
                c.solution.evals_∂ẋ_∂p += 1
            end
            if outputs
                jvp!(c, :∂y_∂p, y_refs, p_refs, p_val, seed; accu = cache.∂y)
                c.solution.evals_∂y_∂p += 1
            end
            if eventIndicators
                jvp!(c, :∂e_∂p, (:indicators, ec_idcs), p_refs, p_val, seed; accu = cache.∂e)
                c.solution.evals_∂e_∂p += 1
            end
        end

        # Time tangent contribution
        if has_t && times && c.fmu.executionConfig.eval_t_gradients
            t_seed = ForwardDiff.partials(t_dual, i)
            if t_seed != 0.0
                if derivatives
                    gvp!(c, :∂ẋ_∂t, dx_refs, :time, t_val, t_seed; accu = cache.∂dx)
                    c.solution.evals_∂ẋ_∂t += 1
                end
                if outputs
                    gvp!(c, :∂y_∂t, y_refs, :time, t_val, t_seed; accu = cache.∂y)
                    c.solution.evals_∂y_∂t += 1
                end
                if eventIndicators
                    gvp!(c, :∂e_∂t, (:indicators, ec_idcs), :time, t_val, t_seed; accu = cache.∂e)
                    c.solution.evals_∂e_∂t += 1
                end
            end
        end

        # Store this partial direction
        @inbounds ∂dx_parts[:, i] .= cache.∂dx
        @inbounds ∂y_parts[:, i] .= cache.∂y
        @inbounds ∂e_parts[:, i] .= cache.∂e
    end

    # --- Assemble Dual output ---
    D = ForwardDiff.Dual{TAG,Float64,N}
    dx_p = cache.dx_primal
    y_p = cache.y_primal
    ec_p = cache.ec_primal

    out_dx = Vector{D}(undef, ndx)
    @inbounds for j in 1:ndx
        out_dx[j] = ForwardDiff.Dual{TAG}(dx_p[j],
            ForwardDiff.Partials(ntuple(i -> ∂dx_parts[j, i], Val(N))))
    end

    out_y = Vector{D}(undef, ny)
    @inbounds for j in 1:ny
        out_y[j] = ForwardDiff.Dual{TAG}(y_p[j],
            ForwardDiff.Partials(ntuple(i -> ∂y_parts[j, i], Val(N))))
    end

    out_ec = Vector{D}(undef, ne)
    @inbounds for j in 1:ne
        out_ec[j] = ForwardDiff.Dual{TAG}(ec_p[j],
            ForwardDiff.Partials(ntuple(i -> ∂e_parts[j, i], Val(N))))
    end

    return FMUEvaluationOutput(out_dx, out_y, out_ec)
end

# ============================================================================
# Entry-point methods for each argument type combination.
# Each extracts primal values and delegates to the shared kernel.
# ============================================================================

# Helper: convert a Real value to Float64
@inline _to_f64(x::Real) = Float64(x)
@inline _to_f64(x::Float64) = x

# Helper: extract Float64 primal value from a ForwardDiff.Dual scalar
@inline _primal_t(t::ForwardDiff.Dual) = Float64(ForwardDiff.value(t))
@inline _primal_t(t::Real) = Float64(t)

# Helper: get cRef and instance
@inline function _get_c_and_cref(cRef)
    cRef_u64 = UInt64(unsense(cRef))
    c = unsafe_pointer_to_objref(Ptr{Nothing}(cRef_u64))
    return c, cRef_u64
end

# Helper: extract Float64 primals from a possibly-Dual vector into a cache buffer
@inline function _get_primals(cache_buf::Vector{Float64}, v::AbstractVector{<:ForwardDiff.Dual})
    _extract_values!(cache_buf, v)
end
@inline _get_primals(::Vector{Float64}, v::AbstractVector{Float64}) = v
@inline _get_primals(::Vector{Float64}, v::AbstractVector{<:Real}) = Float64.(v)

# Macro to reduce boilerplate for ForwardDiff dispatch methods.
# Takes 7 args: (dx_type, y_type, x_type, u_type, p_type, ec_type, t_type)
# Each must be :Dual or :Real, matching the original @ForwardDiff_frule signatures exactly.
macro _fd_eval(dx_type, y_type, x_type, u_type, p_type, ec_type, t_type)
    dx_T = dx_type == :Dual ? :(AbstractVector{<:ForwardDiff.Dual}) : :(AbstractVector{<:Real})
    y_T = y_type == :Dual ? :(AbstractVector{<:ForwardDiff.Dual}) : :(AbstractVector{<:Real})
    x_T = x_type == :Dual ? :(AbstractVector{<:ForwardDiff.Dual{TAG,V,N}}) : :(AbstractVector{<:Real})
    u_T = u_type == :Dual ? :(AbstractVector{<:ForwardDiff.Dual{TAG,V,N}}) : :(AbstractVector{<:Real})
    p_T = p_type == :Dual ? :(AbstractVector{<:ForwardDiff.Dual{TAG,V,N}}) : :(AbstractVector{<:Real})
    ec_T = ec_type == :Dual ? :(AbstractVector{<:ForwardDiff.Dual{TAG,V,N}}) : :(AbstractVector{<:Real})
    t_T = t_type == :Dual ? :(ForwardDiff.Dual{TAG,V,N}) : :Real

    # Primal extraction expressions
    x_primals = x_type == :Dual ? :(_get_primals(fdc.x_val, x)) : :(_as_f64_vec(x))
    u_primals = u_type == :Dual ? :(_get_primals(fdc.u_val, u)) : :(_as_f64_vec(u))
    p_primals = p_type == :Dual ? :(_get_primals(fdc.p_val, p)) : :(_as_f64_vec(p))
    ec_primals = ec_type == :Dual ? :(_get_primals(fdc.ec_val, ec)) : :(_as_f64_vec(ec))
    t_primals = t_type == :Dual ? :(_primal_t(t)) : :(_to_f64(t))

    # Dual pass-through (nothing if not Dual)
    x_duals = x_type == :Dual ? :x : :nothing
    u_duals = u_type == :Dual ? :u : :nothing
    p_duals = p_type == :Dual ? :p : :nothing
    ec_duals = ec_type == :Dual ? :ec : :nothing
    t_duals = t_type == :Dual ? :t : :nothing

    # For dx/y arrays: if they're Dual, allocate Float64 buffers for primal eval
    dx_is_dual = (dx_type == :Dual)
    y_is_dual = (y_type == :Dual)
    dx_expr = dx_is_dual ? :(zeros(Float64, length(dx))) : :dx
    y_expr = y_is_dual ? :(zeros(Float64, length(y))) : :y

    quote
        function FMIBase.eval!(
            cRef::UInt64,
            dx::$dx_T,
            dx_refs::AbstractVector{<:fmiValueReference},
            y::$y_T,
            y_refs::AbstractVector{<:fmiValueReference},
            x::$x_T,
            u::$u_T,
            u_refs::AbstractVector{<:fmiValueReference},
            p::$p_T,
            p_refs::AbstractVector{<:fmiValueReference},
            ec::$ec_T,
            ec_idcs::AbstractVector{<:fmiValueReference},
            t::$t_T,
            x_d::AbstractVector{<:Real},
        ) where {TAG,V,N}
            c, cRef_u64 = _get_c_and_cref(cRef)
            fdc = _get_fd_cache(c)

            dx_refs_clean = _as_u32_vec(dx_refs)
            y_refs_clean = _as_u32_vec(y_refs)
            u_refs_clean = _as_u32_vec(u_refs)
            p_refs_clean = _as_u32_vec(p_refs)
            ec_idcs_clean = _as_u32_vec(ec_idcs)
            x_d_clean = _as_f64_vec(unsense(x_d))

            x_val = $x_primals
            u_val = $u_primals
            p_val = $p_primals
            ec_val = $ec_primals
            t_val = $t_primals

            dx_buf = $dx_expr
            y_buf = $y_expr

            return _eval_forwarddiff!(
                c, cRef_u64,
                dx_buf, dx_refs_clean, y_buf, y_refs_clean,
                x_val, u_val, u_refs_clean, p_val, p_refs_clean,
                ec_val, ec_idcs_clean, t_val, x_d_clean,
                $x_duals, $u_duals, $p_duals, $ec_duals, $t_duals,
                Val{TAG}(), Val{N}(),
            )
        end
    end |> esc
end

# All 17 combinations matching the original @ForwardDiff_frule registrations exactly.
# Format: @_fd_eval dx y x u p ec t

# 1. dx, y, x, u, p, ec, t (all Dual)
@_fd_eval Dual Dual Dual Dual Dual Dual Dual

# 2. dx, y, x, u, t (dx,y,x,u Dual, p,ec Real, t Dual)
@_fd_eval Dual Dual Dual Dual Real Real Dual

# 3. x, u (dx,y Real)
@_fd_eval Real Real Dual Dual Real Real Real

# 4. x, u, t (dx,y Real)
@_fd_eval Real Real Dual Dual Real Real Dual

# 5. x, p
@_fd_eval Real Real Dual Real Dual Real Real

# 6. t only
@_fd_eval Real Real Real Real Real Real Dual

# 7. x only
@_fd_eval Real Real Dual Real Real Real Real

# 8. u only
@_fd_eval Real Real Real Dual Real Real Real

# 9. p only
@_fd_eval Real Real Real Real Dual Real Real

# 10. ec only
@_fd_eval Real Real Real Real Real Dual Real

# 11. x, t
@_fd_eval Real Real Dual Real Real Real Dual

# 12. x, ec, t
@_fd_eval Real Real Dual Real Real Dual Dual

# 13. ec, t
@_fd_eval Real Real Real Real Real Dual Dual

# 14. x, ec
@_fd_eval Real Real Dual Real Real Dual Real

# 15. x, p, t
@_fd_eval Real Real Dual Real Dual Real Dual

# 16. x, p, ec, t
@_fd_eval Real Real Dual Real Dual Dual Dual

# 17. x, p, ec
@_fd_eval Real Real Dual Real Dual Dual Real
