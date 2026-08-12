module MTKOrdinaryDiffEqBDFExt

using ModelingToolkit
using ModelingToolkitBase
using OrdinaryDiffEqBDF
using OrdinaryDiffEqBDF: FBDF, DFBDF
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    probs = ModelingToolkitBase._solver_precompile_problems()
    ode_algs = Any[FBDF()]
    dae_algs = Any[DFBDF()]
    if isdefined(OrdinaryDiffEqBDF, :NordsieckBDF)
        push!(ode_algs, OrdinaryDiffEqBDF.NordsieckBDF())
        push!(dae_algs, OrdinaryDiffEqBDF.DNordsieckBDF())
    end
    @compile_workload begin
        for alg in ode_algs
            ModelingToolkitBase._solve_precompile(probs.ode, alg)
            ModelingToolkitBase._solve_precompile(probs.mmdae, alg)
        end
        for alg in dae_algs
            ModelingToolkitBase._solve_precompile(probs.dae, alg)
        end
    end
end

end
