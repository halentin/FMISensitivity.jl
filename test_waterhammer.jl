using FMIImport, FMIBase

import FMISensitivity.ForwardDiff 
import FMISensitivity.Zygote
import FMISensitivity.ReverseDiff
# import FMISensitivity.FiniteDiff

using FMISensitivity.FMIBase
using FMISensitivity.FMIBase.FMICore

using FMISensitivity.FMIBase: getContinuousStates, getReal, getRealType, getEventIndicators, getDirectionalDerivative

fmu = loadFMU("Waterhammer_massFlowPulse_10V_2.fmu")
md = fmu.modelDescription
c, _ = FMIImport.prepareSolveFMU(fmu, nothing, :ME; loggingOn=true)
# mVs= modelVariablesForValueReference.(Ref(md), md.stateValueReferences)
# names = [v[1].name for v in mVs]
include("dependency_matrix.jl")
dep_mtx = DependencyMatrix(md)
ddx_dx_template = dep_mtx[md.derivativeValueReferences, md.stateValueReferences]
## Only keep "dependent"-Dependencies"
# ddx_dx_template[x->(x>4)]
# ddx_dx_template[findall(x->(x<5), ddx_dx_template)] .=0
# dropzeros(ddx_dx_template)
using SparseDiffTools
colorvec = matrix_colors(ddx_dx_template)
maximum(colorvec)

fmu.executionConfig.eval_t_gradients = true

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

reset! = function(c::FMUInstance)
    c.solution.evals_∂ẋ_∂x = 0
    c.solution.evals_∂ẋ_∂u = 0
    c.solution.evals_∂ẋ_∂p = 0
    c.solution.evals_∂ẋ_∂t = 0

    c.solution.evals_∂y_∂x = 0
    c.solution.evals_∂y_∂u = 0
    c.solution.evals_∂y_∂p = 0
    c.solution.evals_∂y_∂t = 0

    c.solution.evals_∂e_∂x = 0
    c.solution.evals_∂e_∂u = 0
    c.solution.evals_∂e_∂p = 0
    c.solution.evals_∂e_∂t = 0

    # @test length(dx) == length(fmu.modelDescription.derivativeValueReferences)
    # @test length(y) == length(fmu.modelDescription.outputValueReferences)
    # @test length(e) == fmu.modelDescription.numberOfEventIndicators
end

ydx = fmu(;x=x)
ydx = fmu(;x=x, dx_refs=:all)
ydx = fmu(;x=x, u=u, u_refs=u_refs, dx_refs=:all)
ydx = fmu(;x=x, u=u, u_refs=u_refs, y=y, y_refs=y_refs)
ydx = fmu(;dx_refs=:all, x=x, u=u, u_refs=u_refs, y=y, y_refs=y_refs)
ydx = fmu(;x=x, u=u, u_refs=u_refs, y=y, y_refs=y_refs, dx=dx, dx_refs=:all, p=p, p_refs=p_refs)

fmu.executionConfig.JVPBuiltInDerivatives = true
_f = _x -> fmu(;x=_x, dx_refs=:all)
_f(x)
j_fwd = ForwardDiff.jacobian(_f, x)
j_rwd = ReverseDiff.jacobian(_f, x)
j_rwd == zeros(43,43)