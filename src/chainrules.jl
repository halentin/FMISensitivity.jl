#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

function ChainRulesCore.frule(
    Δtuple,
    ::typeof(FMIBase.eval!),
    cRef,
    dx,
    dx_refs,
    y,
    y_refs,
    x,
    u,
    u_refs,
    p,
    p_refs,
    ec,
    ec_idcs,
    t,
    x_d,
)

    Δself,
    ΔcRef,
    Δdx,
    Δdx_refs,
    Δy,
    Δy_refs,
    Δx,
    Δu,
    Δu_refs,
    Δp,
    Δp_refs,
    Δec,
    Δec_idcs,
    Δt,
    Δx_d = Δtuple # undual ?

    @debug "frule start"

    ### ToDo: Somehow, ForwardDiff enters with all types beeing Float64, this needs to be corrected.

    cRef = unsense(cRef) # undual(cRef)
    if typeof(cRef) != UInt64
        cRef = UInt64(cRef)
    end
    c = unsafe_pointer_to_objref(Ptr{Nothing}(cRef))

    # ToDo: is this necessary?
    # t = undual(t)
    # u = undual(u)
    # x = undual(x)
    # p = undual(p)

    dx_refs = unsense(dx_refs)
    dx_refs = convert(Array{UInt32,1}, dx_refs)
    if length(dx_refs) == 0 &&
       length(dx) == length(c.fmu.modelDescription.derivativeValueReferences) # all derivatives, please!
        dx_refs = c.fmu.modelDescription.derivativeValueReferences
    end

    # [Note] `unsense` is necessary for AD over AD
    y_refs = unsense(y_refs)
    u_refs = unsense(u_refs)
    p_refs = unsense(p_refs)
    ec_idcs = unsense(ec_idcs)

    y_refs = convert(Array{UInt32,1}, y_refs)
    u_refs = convert(Array{UInt32,1}, u_refs)
    p_refs = convert(Array{UInt32,1}, p_refs)
    ec_idcs = convert(Array{UInt32,1}, ec_idcs)

    ###

    outputs = (length(y_refs) > 0)
    inputs = (length(u_refs) > 0)
    derivatives = (length(dx) > 0)
    times = FMIBase.isSetReal(c.fmu, t)
    states = (length(x) > 0)
    parameters = (length(p_refs) > 0)
    eventIndicators = (length(ec_idcs) > 0)

    Ω = FMIBase.eval!(
        cRef,
        dx,
        dx_refs,
        y,
        y_refs,
        x,
        u,
        u_refs,
        p,
        p_refs,
        ec,
        ec_idcs,
        t,
        x_d,
    )

    # time, states and inputs where already set in `eval!`, no need to repeat it here

    ∂y = zeros(length(y))
    ∂dx = zeros(length(dx))
    ∂e = zeros(length(ec))

    if Δx != NoTangent() && length(Δx) > 0

        if states
            if derivatives
                jvp!(
                    c,
                    :∂ẋ_∂x,
                    dx_refs,
                    c.fmu.modelDescription.stateValueReferences,
                    x,
                    Δx;
                    accu = ∂dx,
                )
                c.solution.evals_∂ẋ_∂x += 1
            end

            if outputs
                jvp!(
                    c,
                    :∂y_∂x,
                    y_refs,
                    c.fmu.modelDescription.stateValueReferences,
                    x,
                    Δx;
                    accu = ∂y,
                )
                c.solution.evals_∂y_∂x += 1
            end

            if eventIndicators
                jvp!(
                    c,
                    :∂e_∂x,
                    (:indicators, ec_idcs),
                    c.fmu.modelDescription.stateValueReferences,
                    x,
                    Δx;
                    accu = ∂e,
                )
                c.solution.evals_∂e_∂x += 1
            end
        end
    end

    if Δu != NoTangent() && length(Δu) > 0

        if inputs
            if derivatives
                jvp!(c, :∂ẋ_∂u, dx_refs, u_refs, u, Δu; accu = ∂dx)
                c.solution.evals_∂ẋ_∂u += 1
            end

            if outputs
                jvp!(c, :∂y_∂u, y_refs, u_refs, u, Δu; accu = ∂y)
                c.solution.evals_∂y_∂u += 1
            end

            if eventIndicators
                jvp!(c, :∂e_∂u, (:indicators, ec_idcs), u_refs, u, Δu; accu = ∂e)
                c.solution.evals_∂e_∂u += 1
            end
        end
    end

    if Δp != NoTangent() && length(Δp) > 0

        if parameters
            if derivatives
                jvp!(c, :∂ẋ_∂p, dx_refs, p_refs, p, Δp; accu = ∂dx)
                c.solution.evals_∂ẋ_∂p += 1
            end

            if outputs
                jvp!(c, :∂y_∂p, y_refs, p_refs, p, Δp; accu = ∂y)
                c.solution.evals_∂y_∂p += 1
            end

            if eventIndicators
                jvp!(c, :∂e_∂p, (:indicators, ec_idcs), p_refs, p, Δp; accu = ∂e)
                c.solution.evals_∂e_∂p += 1
            end
        end
    end

    if Δt != NoTangent() && c.fmu.executionConfig.eval_t_gradients

        if times
            if derivatives
                gvp!(c, :∂ẋ_∂t, dx_refs, :time, t, Δt; accu = ∂dx)
                c.solution.evals_∂ẋ_∂t += 1
            end

            if outputs
                gvp!(c, :∂y_∂t, y_refs, :time, t, Δt; accu = ∂y)
                c.solution.evals_∂y_∂t += 1
            end

            if eventIndicators
                gvp!(c, :∂e_∂t, (:indicators, ec_idcs), :time, t, Δt; accu = ∂e)
                c.solution.evals_∂e_∂t += 1
            end
        end
    end

    @debug "frule end:   ∂y=$(∂y)   ∂dx=$(∂dx)   ∂e=$(∂e)"

    # [Note] Type Real is required for AD over AD
    ∂Ω = FMUEvaluationOutput{Real}() # Float64
    ∂Ω.dx = ∂dx
    ∂Ω.y = ∂y
    ∂Ω.ec = ∂e

    return Ω, ∂Ω
end

function ChainRulesCore.rrule(
    ::typeof(FMIBase.eval!),
    cRef,
    dx,
    dx_refs,
    y,
    y_refs,
    x,
    u,
    u_refs,
    p,
    p_refs,
    ec,
    ec_idcs,
    t,
    x_d,
)

    @assert !isa(cRef, FMUInstance) "Wrong dispatched!"

    @debug "rrule start: $((cRef, dx, dx_refs, y, y_refs, x, u, u_refs, p, p_refs, ec, ec_idcs, t, x_d))"

    c = unsafe_pointer_to_objref(Ptr{Nothing}(cRef))

    y_len = (isnothing(y_refs) ? 0 : length(y_refs))
    dx_len = (isnothing(dx) ? 0 : length(dx))

    _outputs = (length(y_refs) > 0)
    _derivatives = (length(dx) > 0)
    _eventIndicators = (length(ec) > 0)
    states = (length(x) > 0)
    inputs = (length(u_refs) > 0)
    times = FMIBase.isSetReal(c.fmu, t)
    parameters = (length(p_refs) > 0)

    @assert !issense(x_d) "discrete state sensitive!"

    # two strategies for `snapshotEveryStep`:
    # (false) use the closest snapshot, change values to the current state etc. -> might be difficult with nasty algebraic loops!
    # (true) make snapshots for every time step (more secure, more memory)
    pullback_snapshot = nothing
    Ω = nothing

    if c.fmu.executionConfig.snapshot_every_step

        Ω = FMIBase.eval!(
            cRef,
            dx,
            dx_refs,
            y,
            y_refs,
            x,
            u,
            u_refs,
            p,
            p_refs,
            ec,
            ec_idcs,
            t,
            x_d,
        )

        # [Todo] this is wrong, discrete state may not match, bc rrule could be called after event handling for
        # before a state before the event!
        pullback_snapshot = snapshot!(c)

    else
        Ω = FMIBase.eval!(
            cRef,
            dx,
            dx_refs,
            y,
            y_refs,
            x,
            u,
            u_refs,
            p,
            p_refs,
            ec,
            ec_idcs,
            t,
            x_d,
        )

        # [ToDo] maybe the arrays change between pullback creation and use! check this!
        t = copy(t) # is scalar, but could be AD-primitive.
        x = copy(x)
        x_d = copy(x_d)
        p = copy(p)
        u = copy(u)
    end

    ##############

    if dx_len > 0 && length(dx_refs) == 0 # all derivatives, please!
        dx_refs = c.fmu.modelDescription.derivativeValueReferences
    end
    x_refs = c.fmu.modelDescription.stateValueReferences

    function eval_pullback(r̄)

        @debug "eval pullback start"

        d̄x = @view(r̄[1:dx_len])
        ȳ = @view(r̄[(dx_len+1):(dx_len+y_len)])
        ēc = @view(r̄[(dx_len+y_len+1):end])

        outputs = _outputs && !isZeroTangent(ȳ)
        derivatives = _derivatives && !isZeroTangent(d̄x)
        eventIndicators = _eventIndicators && !isZeroTangent(ēc)

        if !isa(ȳ, AbstractArray)
            ȳ = collect(ȳ)
        end

        if !isa(d̄x, AbstractArray)
            d̄x = collect(d̄x)
        end

        if !isa(ēc, AbstractArray)
            ēc = collect(ēc)
        end

        # here, we need to set the state/time/etc. to fit the instance the pullback was created!
        if c.fmu.executionConfig.snapshot_every_step
            apply!(c, pullback_snapshot)
        end

        # light weight call to eval!
        FMIBase.eval_set!(c, x, u, u_refs, p, p_refs, t, x_d)

        x̄ = zeros(length(x)) #ZeroTangent()
        t̄ = zeros(1) #ZeroTangent()
        ū = zeros(length(u)) #ZeroTangent()
        p̄ = zeros(length(p)) #ZeroTangent()
        x̄_d = zeros(length(x_d)) # ZeroTangent()

        if derivatives
            if states
                vjp!(c, :∂ẋ_∂x, dx_refs, x_refs, x, d̄x; accu = x̄)
                c.solution.evals_∂ẋ_∂x += 1
            end

            if inputs
                vjp!(c, :∂ẋ_∂u, dx_refs, u_refs, u, d̄x; accu = ū)
                c.solution.evals_∂ẋ_∂u += 1
            end

            if parameters
                vjp!(c, :∂ẋ_∂p, dx_refs, p_refs, p, d̄x; accu = p̄)
                c.solution.evals_∂ẋ_∂p += 1
            end

            if times && c.fmu.executionConfig.eval_t_gradients
                vgp!(c, :∂ẋ_∂t, dx_refs, :time, t, d̄x; accu = t̄)
                c.solution.evals_∂ẋ_∂t += 1
            end
        end

        if outputs
            if states
                vjp!(c, :∂y_∂x, y_refs, x_refs, x, ȳ; accu = x̄)
                c.solution.evals_∂y_∂x += 1
            end

            if inputs
                vjp!(c, :∂y_∂u, y_refs, u_refs, u, ȳ; accu = ū)
                c.solution.evals_∂y_∂u += 1
            end

            if parameters
                vjp!(c, :∂y_∂p, y_refs, p_refs, p, ȳ; accu = p̄)
                c.solution.evals_∂y_∂p += 1
            end

            if times && c.fmu.executionConfig.eval_t_gradients
                vgp!(c, :∂y_∂t, y_refs, :time, t, ȳ; accu = t̄)
                c.solution.evals_∂y_∂t += 1
            end
        end

        if _eventIndicators # ToDo: This should be `eventIndicators` but we get it for every ēc bc. of workaround in `condition!`
            if states
                vjp!(c, :∂e_∂x, (:indicators, ec_idcs), x_refs, x, ēc; accu = x̄)
                c.solution.evals_∂e_∂x += 1
            end

            if inputs
                vjp!(c, :∂e_∂u, (:indicators, ec_idcs), u_refs, u, ēc; accu = ū)
                c.solution.evals_∂e_∂u += 1
            end

            if parameters
                vjp!(c, :∂e_∂p, (:indicators, ec_idcs), p_refs, p, ēc; accu = p̄)
                c.solution.evals_∂e_∂p += 1
            end

            if times && c.fmu.executionConfig.eval_t_gradients
                vgp!(c, :∂e_∂t, (:indicators, ec_idcs), :time, t, ēc; accu = t̄)
                c.solution.evals_∂e_∂t += 1
            end
        end

        # write back
        f̄ = [] # NoTangent()
        c̄Ref = [] # ZeroTangent()
        d̄x_refs = [] # ZeroTangent()
        ȳ_refs = [] # ZeroTangent()
        ēc_idcs = [] # ZeroTangent()
        ū_refs = [] # ZeroTangent()
        p̄_refs = [] # ZeroTangent()

        t̄ = t̄[1]

        @debug "pullback on d̄x, ȳ, ēc = $(d̄x), $(ȳ), $(ēc)\nt= $(t)s\nx=$(x)\nx_d=$(x_d)\ndx=$(dx)\n(x̄=$(x̄), x̄_d=$(x̄_d), ū=$(ū), p̄=$(p̄), t̄=$(t̄))"

        if c.fmu.executionConfig.snapshot_every_step
            freeSnapshot!(pullback_snapshot)
        end

        d̄x = zeros(length(dx)) # ZeroTangent()
        ȳ = zeros(length(y)) # ZeroTangent()
        ēc = zeros(length(ec)) # ZeroTangent() # copy(ec) #

        # [ToDo] This needs to be a tuple... but this prevents pre-allocation...
        return (
            f̄,
            c̄Ref,
            d̄x,
            d̄x_refs,
            ȳ,
            ȳ_refs,
            x̄,
            ū,
            ū_refs,
            p̄,
            p̄_refs,
            ēc,
            ēc_idcs,
            t̄,
            x̄_d,
        )
    end

    @debug "rrule end: $((Ω, eval_pullback))"

    return (Ω, eval_pullback)
end