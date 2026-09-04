"""
    precompile_ode_problem()

Small ODE used by the precompile workloads here and in the OrdinaryDiffEq extensions:
tunable parameters, an observed variable eliminated by `mtkcompile`, and a trivial
initialization problem, so the generated problem has the same types as a typical model.
"""
function precompile_ode_problem()
    t = MTKBase.t_nounits
    D = MTKBase.D_nounits
    @parameters a = 1.0 b = 1.0
    @variables x(t) y(t) z(t)
    sys = mtkcompile(
        System([D(x) ~ a * y, D(y) ~ -b * x + z, z ~ x + y], t; name = :precompile_ode)
    )
    return ODEProblem(sys, [x => 1.0, y => 0.0], (0.0, 1.0))
end

"""
    precompile_dae_problem()

Small mass-matrix DAE used by the precompile workloads: one algebraic unknown fixed by a
nonlinear constraint, so `ODEProblem` builds a `NonlinearProblem` for initialization and
`solve` runs `OverrideInit`.
"""
function precompile_dae_problem()
    t = MTKBase.t_nounits
    D = MTKBase.D_nounits
    @parameters k = 1.0
    @variables x(t) y(t)
    sys = mtkcompile(System([D(x) ~ -k * x + y, 0 ~ y^3 + y - x], t; name = :precompile_dae))
    return ODEProblem(sys, [x => 1.0], (0.0, 1.0); guesses = [y => 0.5])
end

"""
    precompile_scc_dae_problem()

Mass-matrix DAE with two nonlinear algebraic unknowns in a chain, so `ODEProblem` builds a
two-block `SCCNonlinearProblem` for initialization.
"""
function precompile_scc_dae_problem()
    t = MTKBase.t_nounits
    D = MTKBase.D_nounits
    @parameters k = 1.0
    @variables x(t)[1:2] y(t)[1:2]
    eqs = [
        D(x[1]) ~ -k * x[1] + y[1], 0 ~ y[1]^3 + y[1] - x[1],
        D(x[2]) ~ -k * x[2] + y[2], 0 ~ y[2]^3 + y[2] - x[2] - y[1],
    ]
    sys = mtkcompile(System(eqs, t; name = :precompile_scc_dae))
    return ODEProblem(
        sys, [sys.x[1] => 1.0, sys.x[2] => 0.5], (0.0, 1.0);
        guesses = [sys.y[1] => 0.5, sys.y[2] => 0.5]
    )
end

"""
    precompile_nlls_problem()

Mass-matrix DAE whose algebraic unknown is also given an initial value, so the
initialization system is overdetermined and `ODEProblem` builds a
`NonlinearLeastSquaresProblem` for it.
"""
function precompile_nlls_problem()
    t = MTKBase.t_nounits
    D = MTKBase.D_nounits
    @parameters k = 1.0
    @variables x(t) y(t)
    sys = mtkcompile(System([D(x) ~ -k * x + y, 0 ~ y^3 + y - x], t; name = :precompile_nlls))
    return ODEProblem(
        sys, [x => 1.0, y => 0.6823278038280193], (0.0, 1.0);
        warn_initialize_determined = false, fully_determined = false
    )
end

"""
    precompile_implicit_dae_problem()

Implicit `DAEProblem` built the way MethodOfLines does it: a first-order system with
algebraic boundary rows, `complete`d but not `mtkcompile`d, with
`build_initializeprob = false` so the solver's own DAE initialization is used.
"""
function precompile_implicit_dae_problem()
    t = MTKBase.t_nounits
    D = MTKBase.D_nounits
    n = 6
    @parameters α = 1.0
    @variables u(t)[1:n]
    dx = 1 / (n - 1)
    eqs = [
        u[1] ~ 0;
        [D(u[i]) ~ α * (u[i + 1] - 2u[i] + u[i - 1]) / dx^2 for i in 2:(n - 1)];
        u[n] ~ 0
    ]
    sys = complete(System(eqs, t; name = :precompile_implicit_dae))
    u0 = [sin(π * (i - 1) * dx) for i in 1:n]
    du0 = [0.0; [(u0[i + 1] - 2u0[i] + u0[i - 1]) / dx^2 for i in 2:(n - 1)]; 0.0]
    op = [[u[i] => u0[i] for i in 1:n]; [D(u[i]) => du0[i] for i in 1:n]]
    return SciMLBase.DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)
end

PrecompileTools.@compile_workload begin
    t = MTKBase.t_nounits
    D = MTKBase.D_nounits

    function f!(du, u, p)
        du[1] = cos(u[2]) - u[1]
        du[2] = sin(u[1] + u[2]) + u[2]
        du[3] = 2u[4] + u[3] + 1.0
        du[4] = u[5]^2 + u[4]
        du[5] = u[3]^2 + u[5]
        du[6] = u[1] + u[2] + u[3] + u[4] + u[5] + 2.0u[6] + 2.5u[7] + 1.5u[8]
        du[7] = u[1] + u[2] + u[3] + 2.0u[4] + u[5] + 4.0u[6] - 1.5u[7] + 1.5u[8]
        du[8] = u[1] + 2.0u[2] + 3.0u[3] + 5.0u[4] + 6.0u[5] + u[6] - u[7] - u[8]

        du[9] = u[4] + u[5] + u[6] + u[7] + u[8] + 5.3u[9] + 5.8u[10] + 4.8u[11]
        du[10] = u[4] + u[5] + u[6] + 5.3u[7] + u[8] + 7.3u[9] - 4.8u[10] + 4.8u[11]
        du[11] = u[4] + 5.3u[5] + 6.3u[6] + 8.3u[7] + 9.3u[8] + u[9] - u[10] - u[11]

        du[12] = u[7] + u[8] + u[9] + u[10] + u[11] + 8.6u[12] + 8.11u[13] + 7.11u[14]
        du[13] = u[7] + u[8] + u[9] + 8.6u[10] + u[11] + 10.6u[12] - 7.11u[13] + 7.11u[14]
        du[14] = u[7] + 8.6u[8] + 9.6u[9] + 11.6u[10] + 12.6u[11] + u[12] - u[13] - u[14]

        du[15] = u[10] + u[11] + u[12] + u[13] + u[14] + 11.9u[15] + 11.14u[16] + 10.14u[17]
        du[16] = u[10] + u[11] + u[12] + 11.9u[13] + u[14] + 13.9u[15] - 10.14u[16] + 10.14u[17]
        du[17] = u[10] + 11.9u[11] + 12.9u[12] + 14.9u[13] + 15.9u[14] + u[15] - u[16] - u[17]

        du[18] = u[13] + u[14] + u[15] + u[16] + u[17] + 14.12u[18] + 14.17u[19] + 13.17u[20]
        du[19] = u[13] + u[14] + u[15] + 14.12u[16] + u[17] + 16.12u[18] - 13.17u[19] + 13.17u[20]
        du[20] = u[13] + 14.12u[14] + 15.12u[15] + 17.12u[16] + 18.12u[17] + u[18] - u[19] - u[20]
    end
    @variables u[1:20] = rand(20) [irreducible = true]
    eqs = Num[0 for _ in 1:20]
    f!(eqs, u, nothing)
    eqs = 0 .~ eqs
    @mtkcompile model = System(eqs)
    sccprob = SCCNonlinearProblem(model, nothing)
    solve(sccprob, SimpleNonlinearSolve.SimpleNewtonRaphson())

    precompile_ode_problem()
    precompile_dae_problem()
    precompile_scc_dae_problem()
    precompile_nlls_problem()
    precompile_implicit_dae_problem()
end
