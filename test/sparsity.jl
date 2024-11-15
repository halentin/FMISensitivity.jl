using FMIImport, FMIZoo
using FMISensitivity
using Test
import FMISensitivity.ForwardDiff
using FMISensitivity
using BenchmarkTools

# using Logging
# debug_logger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debug_logger)

myFMU = loadFMU("waermestab.fmu")
myFMU.executionConfig.use_sparsity_information = false

c, _ = FMIImport.prepareSolveFMU(myFMU, nothing, :ME)

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
myFMU.executionConfig.sensitivity_strategy = :FiniteDiff
j_fwd = ForwardDiff.jacobian(f, x)

@btime ForwardDiff.jacobian(f, x)
function test()
    invalidate_all()
    ForwardDiff.jacobian(f, x)
end
function test()
    invalidate_all()
    FMISensitivity.validate!(c.∂ẋ_∂x, x)
end

@btime test()
unloadFMU(myFMU)

myFMU = loadFMU("waermestab.fmu")
myFMU.executionConfig.use_sparsity_information = true

c, _ = FMIImport.prepareSolveFMU(myFMU, nothing, :ME)

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
myFMU.executionConfig.sensitivity_strategy = :FiniteDiff
j_fwd = ForwardDiff.jacobian(f, x)


colorvec = c.∂ẋ_∂x.default_coloring_col
sparsity = c.∂ẋ_∂x.sparsity_pattern

map = decompression_map(colorvec, sparsity)

@btime ForwardDiff.jacobian(f, x)
function test()
    invalidate_all()
    ForwardDiff.jacobian(f, x)
end
function test()
    invalidate_all()
    FMISensitivity.validate!(c.∂ẋ_∂x, x)
end
@btime test()
unloadFMU(myFMU)