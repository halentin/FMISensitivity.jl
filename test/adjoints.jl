using FMIImport, FMIZoo
using FMISensitivity
import FMISensitivity.validate!
using Test
import FMISensitivity.ForwardDiff


myFMU = loadFMU("VanDerPol", "ModelicaReferenceFMUs", "0.0.30", "3.0")
c, _ = FMIImport.prepareSolveFMU(myFMU, nothing, :ME)
# myFMU.executionConfig.use_sparsity_information = false
x = [2.0, 0.0]
known_result = [[0.0, -1.0] [1.0, -3.0]] # Analytic Solution for the Jacobian at [2.0, 0.0] with μ=1.0
f = x->myFMU(;x=x, dx_refs=:all)
function invalidate_all()
    for sens in [:∂ẋ_∂x, :∂ẋ_∂u, :∂ẋ_∂p, :∂ẋ_∂t, :∂y_∂x, :∂y_∂u, :∂y_∂p, :∂y_∂t, :∂e_∂x, :∂e_∂u, :∂e_∂p, :∂e_∂t]
        field = getfield(c, sens)
        FMISensitivity.invalidate!(field)
    end
end
function delete_all()
    for sens in [:∂ẋ_∂x, :∂ẋ_∂u, :∂ẋ_∂p, :∂ẋ_∂t, :∂y_∂x, :∂y_∂u, :∂y_∂p, :∂y_∂t, :∂e_∂x, :∂e_∂u, :∂e_∂p, :∂e_∂t]
        setfield!(c, sens, nothing)
    end
end

invalidate_all()
delete_all()
myFMU.executionConfig.sensitivity_strategy = :FMIDirectionalDerivative
j_fwd = ForwardDiff.jacobian(f, x)
c.dependency_matrix.matrix
@test j_fwd == known_result

invalidate_all()
delete_all()

myFMU.executionConfig.sensitivity_strategy = :FMIAdjointDerivative
j_fwd = ForwardDiff.jacobian(f, x)
@test j_fwd == known_result


invalidate_all()
delete_all()
myFMU.executionConfig.sensitivity_strategy = :default
j_fwd = ForwardDiff.jacobian(f, x)
@test j_fwd == known_result


# using BenchmarkTools

# function calc_der(f,x)
#     invalidate_all()
#     ForwardDiff.jacobian(f, x)
# end

# delete_all()
# myFMU.executionConfig.sensitivity_strategy = :FMIDirectionalDerivative
# @btime calc_der($f,$x)
# delete_all()
# myFMU.executionConfig.sensitivity_strategy = :FMIAdjointDerivative
# @btime calc_der($f,$x)
# delete_all()
# myFMU.executionConfig.sensitivity_strategy = :default
# @btime calc_der($f,$x)


# function test_validate(c, x)
#     invalidate_all()
#     FMISensitivity.validate!(c.∂ẋ_∂x, x)
# end
# invalidate_all()
# test_validate()
# @btime test_validate($c, $x)