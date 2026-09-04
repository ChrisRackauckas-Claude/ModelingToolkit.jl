module MTKOrdinaryDiffEqBDFExt

using ModelingToolkit
using OrdinaryDiffEqBDF: FBDF, DFBDF
using DiffEqBase: BrownFullBasicInit
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    odeprob = ModelingToolkit.precompile_ode_problem()
    daeprob = ModelingToolkit.precompile_dae_problem()
    sccprob = ModelingToolkit.precompile_scc_dae_problem()
    nllsprob = ModelingToolkit.precompile_nlls_problem()
    implicitprob = ModelingToolkit.precompile_implicit_dae_problem()
    @compile_workload begin
        solve(odeprob, FBDF())
        solve(daeprob, FBDF())
        solve(sccprob, FBDF())
        solve(nllsprob, FBDF())
        solve(implicitprob, DFBDF(); initializealg = BrownFullBasicInit())
    end
end

end
