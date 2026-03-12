#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

# FiniteDiff Jacobians

abstract type FMUSensitivities end

mutable struct FMUJacobian{C,T,F} <: FMUSensitivities
    valid::Bool
    colored::Bool
    instance::C

    mtx::Matrix{T}
    jvp::Vector{T}
    vjp::Vector{T}

    f_refs::Union{Vector{UInt32},Tuple{Symbol,Vector{UInt32}}}
    x_refs::Union{Vector{UInt32},Symbol}
    f_refs_set::Union{Set,Nothing}

    f::F

    #cache::FiniteDiff.JacobianCache
    #colors::

    validations::Int
    colorings::Int

    seed::Vector{T}   # pre-allocated seed for in-place partial extraction

    function FMUJacobian{T}(
        instance::C,
        f_refs::Union{Vector{UInt32},Tuple{Symbol,Vector{UInt32}}},
        x_refs::Union{Vector{UInt32},Symbol},
    ) where {C,T}

        @assert !isa(f_refs, Tuple) || f_refs[1] == :indicators "`f_refs` is Tuple, it must be `:indicators`"
        @assert !isa(x_refs, Symbol) || x_refs == :time "`x_refs` is Symbol, it must be `:time`"

        f_len = 0
        x_len = 0
        f_refs_set = nothing
        f = nothing

        if isa(f_refs, Tuple)
            f_len = length(f_refs[2]) # number of event indicators to capture
            x_len = length(x_refs)
            f = f_∂e_∂v
        else
            f_len = length(f_refs)
            x_len = length(x_refs)
            f_refs_set = Set(f_refs)
            f = f_∂v_∂v
        end

        F = typeof(f)

        inst = new{C,T,F}()
        inst.f = f
        inst.instance = instance
        inst.f_refs = f_refs
        inst.f_refs_set = f_refs_set
        inst.x_refs = x_refs

        inst.mtx = zeros(T, f_len, x_len)
        inst.jvp = zeros(T, f_len)
        inst.vjp = zeros(T, x_len)

        inst.valid = false
        inst.validations = 0
        inst.colored = false
        inst.colorings = 0

        inst.seed = zeros(T, x_len)

        return inst
    end

end

mutable struct FMUGradient{C,T,F} <: FMUSensitivities
    valid::Bool
    colored::Bool
    instance::C

    vec::Vector{T}
    gvp::Vector{T}
    vgp::Vector{T}

    f_refs::Union{Vector{UInt32},Tuple{Symbol,Vector{UInt32}}}
    x_refs::Union{Vector{UInt32},Symbol}
    f_refs_set::Union{Set,Nothing}

    f::F

    #cache::FiniteDiff.GradientCache
    #colors::

    validations::Int
    colorings::Int

    seed::Vector{T}   # pre-allocated seed for in-place partial extraction

    function FMUGradient{T}(
        instance::C,
        f_refs::Union{Vector{UInt32},Tuple{Symbol,Vector{UInt32}}},
        x_refs::Union{UInt32,Symbol},
    ) where {C,T}

        @assert !isa(f_refs, Tuple) || f_refs[1] == :indicators "`f_refs` is Tuple, it must be `:indicators`"
        @assert !isa(x_refs, Symbol) || x_refs == :time "`x_refs` is Symbol, it must be `:time`"

        f_len = 0
        x_len = 1
        f_refs_set = nothing
        f = nothing

        if isa(f_refs, Tuple)
            f_len = length(f_refs[2])
            f = f_∂e_∂t
        else
            f_len = length(f_refs)
            f_refs_set = Set(f_refs)
            f = f_∂v_∂t
        end

        F = typeof(f)

        inst = new{C,T,F}()
        inst.f = f
        inst.instance = instance
        inst.f_refs = f_refs
        inst.f_refs_set = f_refs_set
        inst.x_refs = x_refs

        inst.vec = zeros(T, f_len)
        inst.gvp = zeros(T, f_len)
        inst.vgp = zeros(T, x_len)

        inst.valid = false
        inst.validations = 0
        inst.colored = false
        inst.colorings = 0

        inst.seed = zeros(T, 1)

        return inst
    end

end

function f_∂v_∂v(jac::FMUJacobian, f, x)
    setReal(jac.instance, jac.x_refs, x; track = false)
    getReal!(jac.instance, jac.f_refs, f)
    return f
end

function f_∂e_∂v(jac::FMUJacobian, f, x)
    symbol, f_refs = jac.f_refs
    @assert symbol == :indicators "Called `f_∂e_∂v` but f_refs is not in event indicator shape."
    setReal(jac.instance, jac.x_refs, x; track = false)
    getEventIndicators!(jac.instance, f, f_refs)
    return f
end

function f_∂e_∂t(jac::FMUGradient, f, x)
    setTime(jac.instance, x; track = false)
    getEventIndicators!(jac.instance, f, jac.f_refs[2])
    return f
end

function f_∂v_∂t(jac::FMUGradient, f, x)
    setTime(jac.instance, x; track = false)
    getReal!(jac.instance, jac.f_refs, f)
    return f
end

function FMIBase.invalidate!(sens::FMUSensitivities)
    sens.valid = false
    return nothing
end

function FMIBase.check_invalidate!(vrs, sens::FMUSensitivities)
    if !sens.valid
        return
    end

    if isnothing(sens.f_refs_set)
        return
    end

    for vr ∈ vrs
        if vr ∈ sens.f_refs_set
            invalidate!(sens)
        end
    end

    return nothing
end

function uncolor!(jac::FMUSensitivities)
    jac.colored = false
    return nothing
end

function onehot!(seed, i::Integer) # [ToDo] this could be solved without allocations
    seed .= 0.0
    seed[i] = 1.0
    return seed
end

function validate!(jac::FMUJacobian, x::AbstractVector)

    rows = length(jac.f_refs)
    cols = length(jac.x_refs)

    # only VR to VR value references can be sampled using built-in functions in FMI
    if !isa(jac.f_refs, Tuple) && !isa(jac.x_refs, Symbol)
        if jac.instance.fmu.executionConfig.sensitivity_strategy ==
           :FMIDirectionalDerivative && providesDirectionalDerivatives(jac.instance.fmu)

            # ToDo: use directional derivatives with sparsitiy information!
            # ToDo: Optimize allocation (onehot)
            # [Note] Jacobian is sampled column by column

            seed = zeros(getRealType(jac.instance), cols)

            for i = 1:cols
                status = getDirectionalDerivative!(
                    jac.instance,
                    jac.f_refs,
                    jac.x_refs,
                    onehot!(seed, i),
                    view(jac.mtx, 1:rows, i),
                )
            end
        elseif jac.instance.fmu.executionConfig.sensitivity_strategy ==
               :FMIAdjointDerivative && providesAdjointDerivatives(jac.instance.fmu)

            # ToDo: use directional derivatives with sparsitiy information!
            # ToDo: Optimize allocation (onehot)
            # [Note] Jacobian is sampled row by row

            seed = zeros(getRealType(jac.instance), rows)

            for i = 1:rows
                getAdjointDerivative!(
                    jac.instance,
                    jac.f_refs,
                    jac.x_refs,
                    onehot!(seed, i),
                    view(jac.mtx, 1:cols, i),
                )
            end
        elseif jac.instance.fmu.executionConfig.sensitivity_strategy == :FiniteDiff

            seed = zeros(getRealType(jac.instance), cols)

            # ToDo: also use FiniteDiff here!
            #finite_diff_jacobian!(jac, x)

            for i = 1:cols
                sampleDirectionalDerivative!(
                    jac.instance,
                    jac.f_refs,
                    jac.x_refs,
                    onehot!(seed, i),
                    view(jac.mtx, 1:rows, i);
                    Δx = jac.instance.fmu.executionConfig.finitediff_absstep,
                )
            end
        else
            @assert false "Unknown sensitivity strategy `$(jac.instance.fmu.executionConfig.sensitivity_strategy)`."
        end
    else
        finite_diff_jacobian!(jac, x)
    end

    jac.validations += 1
    jac.valid = true
    return nothing
end

function finite_diff_jacobian!(jac, x)

    # FMUs remember their state, therefore me need to check the state before sampling ...
    if !isa(jac.x_refs, Symbol)
        x_old = FMIBase.getReal(jac.instance, jac.x_refs)
    end

    # cache = FiniteDiff.JacobianCache(x)
    fdtype = jac.instance.fmu.executionConfig.finitediff_fdtype

    # this is FiniteDiff default behaviour
    relstep = FiniteDiff.default_relstep(fdtype, eltype(x))
    absstep = relstep

    if jac.instance.fmu.executionConfig.finitediff_relstep >= 0.0
        relstep = jac.instance.fmu.executionConfig.finitediff_relstep
    end

    if jac.instance.fmu.executionConfig.finitediff_absstep >= 0.0
        absstep = jac.instance.fmu.executionConfig.finitediff_absstep
    end

    #@info "x: $(x)"
    #@info "size(jac.mtx): $(size(jac.mtx))"

    #jac.mtx = transpose(jac.mtx)

    # ToDo: for setting `fdtype`, a weird error message is generated (looks like a different sampling pattern)
    FiniteDiff.finite_difference_jacobian!(
        jac.mtx,
        (_dx, _x) -> jac.f(jac, _dx, _x),
        x,
        fdtype;
        relstep = relstep,
        absstep = absstep, #
    ) # , cache)

    #jac.mtx = transpose(jac.mtx)

    # ... and set it afterwards
    if !isa(jac.x_refs, Symbol)
        FMIBase.setReal(jac.instance, jac.x_refs, x_old)
    end
    return nothing
end

function finite_diff_gradient!(grad, x)

    # FMUs remember their state, therefore me need to check the state before sampling ...
    if !isa(grad.x_refs, Symbol)
        x_old = FMIBase.getReal(grad.instance, grad.x_refs)
    end

    # cache = FiniteDiff.JacobianCache(x)
    fdtype = grad.instance.fmu.executionConfig.finitediff_fdtype

    # this is FiniteDiff default behaviour
    relstep = FiniteDiff.default_relstep(fdtype, eltype(x))
    absstep = relstep

    if grad.instance.fmu.executionConfig.finitediff_relstep >= 0.0
        relstep = grad.instance.fmu.executionConfig.finitediff_relstep
    end

    if grad.instance.fmu.executionConfig.finitediff_absstep >= 0.0
        absstep = grad.instance.fmu.executionConfig.finitediff_absstep
    end

    # cache = FiniteDiff.GradientCache(x)
    FiniteDiff.finite_difference_gradient!(
        grad.vec,
        (_dx, _x) -> (grad.f(grad, _dx, _x)),
        x,
        fdtype;
        relstep = relstep,
        absstep = absstep,
    ) # , cache)

    # ... and set it afterwards
    if !isa(grad.x_refs, Symbol)
        FMIBase.setReal(grad.instance, grad.x_refs, x_old)
    end
    return nothing
end

function validate!(grad::FMUGradient, x::Real)

    if !isa(grad.f_refs, Tuple) && !isa(grad.x_refs, Symbol)

        if grad.instance.fmu.executionConfig.sensitivity_strategy ==
           :FMIDirectionalDerivative && providesDirectionalDerivatives(grad.instance.fmu)

            # ToDo: use directional derivatives with sparsitiy information!
            getDirectionalDerivative!(
                grad.instance,
                grad.f_refs,
                grad.x_refs,
                ones(length(jac.f_refs)),
                grad.vec,
            )
        elseif grad.instance.fmu.executionConfig.sensitivity_strategy == :FiniteDiff
            finite_diff_gradient!(grad, x)
        else
            @assert false "Unknown sensitivity strategy `$(grad.instance.fmu.executionConfig.sensitivity_strategy)`."
        end
    else
        finite_diff_gradient!(grad, x)
    end

    grad.validations += 1
    grad.valid = true
    return nothing
end

function color!(sens::FMUSensitivities)
    # ToDo
    # colors = SparseDiffTools.matrix_colors(sparsejac)

    sens.colorings += 1
    sens.colored = true
    return nothing
end

function ref_length(ref::AbstractArray)
    return length(ref)
end

function ref_length(ref::Symbol)
    if ref == :time
        return 1
    else
        @assert false "unknwon ref symbol: $(ref)"
    end
end

function ref_length(ref::Tuple)
    @assert length(ref) == 2 "tuple ref length is $(length(ref)) != 2"
    if ref[1] == :indicators
        return length(ref[2])
    else
        @assert false "unknwon tuple ref $(ref)"
    end
end

function update!(jac::FMUJacobian, x)

    if size(jac.mtx) != (ref_length(jac.f_refs), ref_length(jac.x_refs))
        #if length(jac.mtx) != ref_length(jac.f_refs) * ref_length(jac.x_refs) # this is cheaper
        jac.mtx = similar(jac.mtx, ref_length(jac.f_refs), ref_length(jac.x_refs))
        jac.jvp = similar(jac.jvp, ref_length(jac.f_refs))
        jac.vjp = similar(jac.vjp, ref_length(jac.x_refs))

        jac.valid = false
    end

    if !jac.valid
        validate!(jac, x)
    end

    if !jac.colored
        color!(jac)
    end
    return nothing
end

function update!(gra::FMUGradient, x)

    if length(gra.vec) != ref_length(gra.f_refs)
        gra.vec = similar(gra.vec, ref_length(gra.f_refs))
        gra.gvp = similar(gra.gvp, ref_length(gra.f_refs))
        gra.vgp = similar(gra.vgp, ref_length(gra.x_refs))

        gra.valid = false
    end

    if !gra.valid
        validate!(gra, x)
    end

    if !gra.colored
        color!(gra)
    end
    return nothing
end

function jvp!(jac::FMUJacobian, x::AbstractVector, v::AbstractVector; jvp = jac.jvp)
    FMISensitivity.update!(jac, x)
    #return jac.mtx * v
    mul!(jvp, jac.mtx, v)
    return nothing
end

function vjp!(jac::FMUJacobian, x::AbstractVector, v::AbstractVector; vjp = jac.vjp)
    FMISensitivity.update!(jac, x)
    #return jac.mtx' * v
    mul!(vjp, jac.mtx', v)
    return nothing
end

function gvp!(grad::FMUGradient, x, v; gvp = grad.gvp)
    FMISensitivity.update!(grad, x)
    #return grad.vec * v
    mul!(gvp, grad.vec, v)
    return nothing
end

function vgp!(grad::FMUGradient, x, v, vgp = grad.vgp)
    FMISensitivity.update!(grad, x)
    mul!(vgp, grad.vec', v)
    return nothing
end
