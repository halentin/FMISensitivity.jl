#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

import FMIBase: eval!, invalidate!, check_invalidate!
using FMIBase:
    getDirectionalDerivative!, getAdjointDerivative!, sampleDirectionalDerivative!
using FMIBase:
    setContinuousStates,
    setInputs,
    setReal,
    setTime,
    setReal,
    getReal!,
    getEventIndicators!,
    getRealType,
    startSampling,
    stopSampling,
    issense

# in FMI2 and FMI3 we can use fmi2GetDirectionalDerivative for JVP-computations
function jvp!(c::FMUInstance, mtxCache::Symbol, ∂f_refs, ∂x_refs, x, seed; accu = nothing)

    jac = getfield(c, mtxCache)
    if isnothing(jac)
        # [Note] type Real, so AD-primitves can be stored for AD over AD
        # this is necessary for e.g. gradient over implicit solver solutions with autodiff=true
        T = typeof(seed[1])
        jac = FMUJacobian{T}(c, ∂f_refs, ∂x_refs)
        setfield!(c, mtxCache, jac)
    end

    jac.f_refs = ∂f_refs
    jac.x_refs = ∂x_refs

    if c.fmu.executionConfig.JVPBuiltInDerivatives &&
       providesDirectionalDerivatives(c.fmu) &&
       !isa(jac.f_refs, Tuple) &&
       !isa(jac.x_refs, Symbol)
        getDirectionalDerivative!(c, ∂f_refs, ∂x_refs, seed, jac.jvp)
    else
        jvp!(jac, x, seed)
    end

    accu .+= jac.jvp

    return nothing
end

function gvp!(c::FMUInstance, mtxCache::Symbol, ∂f_refs, ∂x_refs, x, seed; accu = nothing)

    grad = getfield(c, mtxCache)
    if isnothing(grad)
        # [Note] type Real, so AD-primitves can be stored for AD over AD
        # this is necessary for e.g. gradient over implicit solver solutions with autodiff=true
        T = typeof(seed[1])
        grad = FMUGradient{T}(c, ∂f_refs, ∂x_refs)
        setfield!(c, mtxCache, grad)
    end

    grad.f_refs = ∂f_refs
    grad.x_refs = ∂x_refs

    if c.fmu.executionConfig.JVPBuiltInDerivatives &&
       providesDirectionalDerivatives(c.fmu) &&
       !isa(grad.f_refs, Tuple) &&
       !isa(grad.x_refs, Symbol)
        getDirectionalDerivative!(c, ∂f_refs, ∂x_refs, [seed], grad.gvp)
    else
        gvp!(grad, x, seed)
    end

    accu .+= grad.gvp

    return nothing
end

# in FMI2 there is no helper for VJP-computations (but in FMI3) ...
function vjp!(c::FMUInstance, mtxCache::Symbol, ∂f_refs, ∂x_refs, x, seed; accu = nothing)

    jac = getfield(c, mtxCache)
    if isnothing(jac)
        # [Note] type Real, so AD-primitves can be stored for AD over AD
        # this is necessary for e.g. gradient over implicit solver solutions with autodiff=true
        T = typeof(seed[1])
        jac = FMUJacobian{T}(c, ∂f_refs, ∂x_refs)
        setfield!(c, mtxCache, jac)
    end

    jac.f_refs = ∂f_refs
    jac.x_refs = ∂x_refs

    if c.fmu.executionConfig.VJPBuiltInDerivatives &&
       providesAdjointDerivatives(c.fmu) &&
       !isa(jac.f_refs, Tuple) &&
       !isa(jac.x_refs, Symbol)
        getAdjointDerivative!(c, ∂f_refs, ∂x_refs, seed, jac.vjp)
    else
        vjp!(jac, x, seed)
    end

    accu .+= jac.vjp

    return nothing
end

function vgp!(c::FMUInstance, mtxCache::Symbol, ∂f_refs, ∂x_refs, x, seed; accu = nothing)

    grad = getfield(c, mtxCache)
    if isnothing(grad)
        # [Note] type Real, so AD-primitves can be stored for AD over AD
        # this is necessary for e.g. gradient over implicit solver solutions with autodiff=true
        T = typeof(seed[1])
        grad = FMUGradient{T}(c, ∂f_refs, ∂x_refs)
        setfield!(c, mtxCache, grad)
    end

    grad.f_refs = ∂f_refs
    grad.x_refs = ∂x_refs

    if c.fmu.executionConfig.VJPBuiltInDerivatives &&
       providesAdjointDerivatives(c.fmu) &&
       !isa(grad.f_refs, Tuple) &&
       !isa(grad.x_refs, Symbol)
        getAdjointDerivative!(c, ∂f_refs, ∂x_refs, [seed], grad.vgp)
    else
        vgp!(grad, x, seed)
    end

    accu .+= grad.vgp

    return nothing
end
