using FMIImport, FMIZoo
using FMISensitivity
using Test
import FMISensitivity.ForwardDiff
using BenchmarkTools

myFMU = loadFMU("waermestab.fmu")
c, _ = FMIImport.prepareSolveFMU(myFMU, nothing, :ME)

md
x = zeros(16)
f = x->myFMU(;x=x, dx_refs=:all)
function invalidate_all()
    for sens in [:∂ẋ_∂x, :∂ẋ_∂u, :∂ẋ_∂p, :∂ẋ_∂t, :∂y_∂x, :∂y_∂u, :∂y_∂p, :∂y_∂t, :∂e_∂x, :∂e_∂u, :∂e_∂p, :∂e_∂t]
        field = getfield(c, sens)
        FMISensitivity.invalidate!(field)
    end
end
invalidate_all()
myFMU.executionConfig.sensitivity_strategy = :FMIDirectionalDerivative
j_fwd = ForwardDiff.jacobian(f, x)

@btime ForwardDiff.jacobian(f, x)
function test()
    invalidate_all()
    ForwardDiff.jacobian(f, x)
end
@btime test()
