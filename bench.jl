import Pkg
Pkg.activate(joinpath(@__DIR__, "test"))

using FMISensitivity
using FMISensitivity: ForwardDiff, ReverseDiff, FiniteDiff
using FMISensitivity.FMIBase
using FMISensitivity.FMIBase.FMICore
using FMISensitivity.FMIBase: getContinuousStates, getReal, getRealType, getEventIndicators

using FMIImport
using FMIZoo
using Printf

function bench(name, f; N=100)
    f() # extra warmup
    times = Float64[]
    allocs_list = Int[]
    for _ in 1:N
        stats = @timed f()
        push!(times, stats.time)
        push!(allocs_list, Int(stats.bytes))
    end
    sort!(times)
    sort!(allocs_list)
    med_time = times[N÷2]
    med_allocs = allocs_list[N÷2]
    @printf("  %-45s %8.1f μs  %8d bytes\n", name, med_time * 1e6, med_allocs)
end

function run_model_bench(name, tool, year, fmi_ver)
    fmu = loadFMU(name, tool, year, fmi_ver)
    fmu.executionConfig.eval_t_gradients = true
    c, _ = FMIImport.prepareSolveFMU(fmu, nothing, :ME; loggingOn=false)

    x = getContinuousStates(c)
    dx_refs = c.fmu.modelDescription.derivativeValueReferences
    dx = getReal(c, dx_refs)
    u_refs = fmu.modelDescription.inputValueReferences
    u = ones(getRealType(fmu), max(1, length(u_refs)))
    if isempty(u_refs)
        u_refs = fmu.modelDescription.stateValueReferences[1:1]
    end
    y_refs = fmu.modelDescription.outputValueReferences
    p_refs = fmu.modelDescription.parameterValueReferences
    p = getReal(c, p_refs)
    ec_idcs = collect(UInt32(i) for i in 1:fmu.modelDescription.numberOfEventIndicators)
    t = 0.0

    nx = length(x)
    ny = length(y_refs)
    nec = length(ec_idcs)
    np = length(p)

    println("\n=== $name (states=$nx, outputs=$ny, events=$nec, params=$np) ===\n")

    # Test functions
    f_jac_x = _x -> fmu(; x=_x, dx_refs=:all, y_refs=y_refs).buffer
    f_grad_t = _t -> fmu(; t=_t, dx_refs=:all)

    has_ec = nec > 0
    f_jac_xec = has_ec ? (_x -> fmu(; x=_x, dx_refs=:all, y_refs=y_refs, ec_idcs=ec_idcs).buffer) : nothing
    f_jac_p = np > 0 ? (_p -> fmu(; p=_p, p_refs=p_refs, dx_refs=:all, y_refs=y_refs).buffer) : nothing

    # Warmup
    for _ in 1:2
        f_jac_x(x); ForwardDiff.jacobian(f_jac_x, x); ReverseDiff.jacobian(f_jac_x, x)
        if has_ec
            f_jac_xec(x); ForwardDiff.jacobian(f_jac_xec, x); ReverseDiff.jacobian(f_jac_xec, x)
        end
        if f_jac_p !== nothing
            f_jac_p(p); ForwardDiff.jacobian(f_jac_p, p); ReverseDiff.jacobian(f_jac_p, p)
        end
        f_grad_t(t); ForwardDiff.derivative(f_grad_t, t)
    end

    println("ForwardDiff:")
    bench("jacobian ∂(dx,y)/∂x     [$(nx)→$(nx+ny)]", () -> ForwardDiff.jacobian(f_jac_x, x))
    if has_ec
        bench("jacobian ∂(dx,y,ec)/∂x [$(nx)→$(nx+ny+nec)]", () -> ForwardDiff.jacobian(f_jac_xec, x))
    end
    if f_jac_p !== nothing
        bench("jacobian ∂(dx,y)/∂p     [$(np)→$(nx+ny)]", () -> ForwardDiff.jacobian(f_jac_p, p))
    end
    bench("derivative ∂dx/∂t       [1→$(nx)]", () -> ForwardDiff.derivative(f_grad_t, t))

    println("\nReverseDiff:")
    bench("jacobian ∂(dx,y)/∂x     [$(nx)→$(nx+ny)]", () -> ReverseDiff.jacobian(f_jac_x, x))
    if has_ec
        bench("jacobian ∂(dx,y,ec)/∂x [$(nx)→$(nx+ny+nec)]", () -> ReverseDiff.jacobian(f_jac_xec, x))
    end
    if f_jac_p !== nothing
        bench("jacobian ∂(dx,y)/∂p     [$(np)→$(nx+ny)]", () -> ReverseDiff.jacobian(f_jac_p, p))
    end

    unloadFMU(fmu)
end

# Small model
run_model_bench("SpringFrictionPendulumExtForce1D", "Dymola", "2023x", "2.0")

# Model with many params (12 params → wider Jacobian)
run_model_bench("SpringFrictionPendulum1D", "Dymola", "2023x", "2.0")

# Model with many events (34 event indicators)
run_model_bench("SpringTimeFrictionPendulum1D", "Dymola", "2023x", "2.0")
