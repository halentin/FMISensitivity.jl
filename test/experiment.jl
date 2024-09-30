import FMISensitivity.ForwardDiff 
import FMISensitivity.Zygote
import FMISensitivity.ReverseDiff
# import FMISensitivity.FiniteDiff
using FMISensitivity
using FMISensitivity.FMIBase
using FMISensitivity.FMIBase.FMICore

using FMISensitivity.FMIBase: getContinuousStates, getReal, getRealType, getEventIndicators, getDirectionalDerivative


ENV["EXPORTINGTOOL"] = "Dymola"
ENV["EXPORTINGVERSION"] = "2023x"
ENV["FMIVERSION"] = 3.0
ENV["FMUSTRUCT"] = "FMUCOMPONENT"

c, fmu = getFMUStruct("SpringFrictionPendulumExtForce1D", :ME)

fmu.modelDescription

x_refs = fmu.modelDescription.stateValueReferences
x = getContinuousStates(c) 
dx_refs = c.fmu.modelDescription.derivativeValueReferences
dx = getReal(c, dx_refs)
u_refs = fmu.modelDescription.inputValueReferences
u = zeros(getRealType(fmu), length(u_refs))
y_refs = fmu.modelDescription.outputValueReferences
y = getReal(c, y_refs)
p_refs = fmu.modelDescription.parameterValueReferences
p = getReal(c, p_refs)
e = getEventIndicators(c)
t = 0.0

import FMISensitivity.FiniteDiff
_f = _x -> fmu(;x=_x, dx_refs=:all)
_f = _x -> convert(Vector{Real}, fmu(;x=_x, dx_refs=:all))
_f! = (fx, _x) -> fx = fmu(;x=_x, dx_refs=:all)
jac_cache = FiniteDiff.JacobianCache(x)
output = zeros(2,2)
jacs = FiniteDiff.finite_difference_jacobian!(output, _f!, x, jac_cache)
jacs = FiniteDiff.finite_difference_jacobian(_f, x, jac_cache)
j_fwd = ForwardDiff.jacobian(_f, x)

f_test = x -> [x[1]^2+x[2], x[2]^2]
FiniteDiff.finite_difference_jacobian(f_test, [1.0, 2.0])
f_test([2.0 1.0])

a = _f(x)
using FiniteDifferences
_f(x)
method = central_fdm(2,1)
jacobian(method,_f,x)[1]
a = _f(x)

_f = _p -> convert(Vector{Real},fmu(;p=_p, p_refs=p_refs, y_refs=y_refs))
_f(p)
j_fwd = ForwardDiff.jacobian(_f, p)
j_rwd = ReverseDiff.jacobian(_f, p)

j_fd = jacobian(method, _f, p)[1]

_f = _u -> fmu(; u=_u, u_refs=u_refs, dx_refs=:all)
_f(u)
j_fd = jacobian(method, _f, u)[1]


_f = _p -> fmu(;p=_p, p_refs=p_refs, dx_refs=:all)
@test _f(p-ones(12)) == _f(p+ones(12))

fmu.executionConfig.sensitivity_strategy = :FMIDirectionalDerivative
f = x->fmu(;x=x, dx_refs=:all)
j_fwd = ForwardDiff.jacobian(f, x)

invalidate_all()
fmu.executionConfig.sensitivity_strategy = :FiniteDiff
f = _p->fmu(;p=_p, p_refs=p_refs, dx_refs=:all)
j_fwd = ForwardDiff.jacobian(f, p)

invalidate_all()
fmu.executionConfig.sensitivity_strategy = :FMIDirectionalDerivative
f = _p->fmu(;p=_p, p_refs=p_refs, dx_refs=:all)
f(p)
j_fwd = ForwardDiff.jacobian(f, p)
∂ẋ_∂p = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0375 0.0 0.0 10.0 0.6 10.0 0.0 0.0 0.0 5.0 -5.25 0.0]
@test isapprox(j_fwd, ∂ẋ_∂p, atol = 0.001)
function invalidate_all()
    for sens in [:∂ẋ_∂x, :∂ẋ_∂u, :∂ẋ_∂p, :∂ẋ_∂t, :∂y_∂x, :∂y_∂u, :∂y_∂p, :∂y_∂t, :∂e_∂x, :∂e_∂u, :∂e_∂p, :∂e_∂t]
        field = getfield(c, sens)
        FMISensitivity.invalidate!(field)
    end
end
using BenchmarkTools

result = zeros(2)
using FMIBase: getDirectionalDerivative!

function onehot(c::FMUInstance, len::Integer, i::Integer) # [ToDo] this could be solved without allocations
    ret = zeros(getRealType(c), len)
    ret[i] = 1.0
    return ret 
end
using LinearAlgebra
Id = Matrix{Float64}(I, 2,2)
Id[:,1]
view(Id, 1:2,2)
a = onehot(c, 2, 1)
ve = [0.0, 1.0, 0.0]
p_ve = pointer(ve)
@btime sum(view(ve, 1:2))
a = [1.0]

function test_1(c, dx_refs, x_refs, result_mtx)
    rows = length(dx_refs)
    cols = length(x_refs)
    for i in 1:rows
        getDirectionalDerivative!(c, dx_refs, x_refs, onehot(c,cols,i), view(result_mtx, :, i))
    end
end

function test_2(c, dx_refs, x_refs, result_mtx)
    rows = length(dx_refs)
    cols = length(x_refs)
    cache = zeros(cols)
    for i in 1:rows
        getDirectionalDerivative!(c, dx_refs, x_refs, onehot(c,cols,i), cache)
        result_mtx[:,i] = cache
    end
end

function test_3(c, dx_refs, x_refs, result_mtx, perturbation)
    rows = length(dx_refs)
    cols = length(x_refs)
    for i in 1:rows
        v_xrefs = view(x_refs,i)
        v_results = view(result_mtx, :, i)
        getDirectionalDerivative!(c, dx_refs, v_xrefs, perturbation, v_results)
    end
end

result = zeros(2,2)

@btime test_1(c, dx_refs, x_refs, result)
@btime test_2(c, dx_refs, x_refs, result)
perturbation = [1.0]
@btime test_3(c, dx_refs, x_refs, result, perturbation)
@btime getDirectionalDerivative!(c, dx_refs, x_refs, view(ve, 1:2), result)
@btime getDirectionalDerivative!(c, dx_refs, view(x_refs, 1:1), a, result)
@btime getDirectionalDerivative!(c, dx_refs, x_refs, p_ve, result)


result