module MTKOrdinaryDiffEqDefaultExt

using ModelingToolkit
using ModelingToolkitBase
using OrdinaryDiffEqDefault: OrdinaryDiffEqDefault
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    probs = ModelingToolkitBase._solver_precompile_problems()
    @compile_workload begin
        solve(probs.ode; save_everystep = false)
        solve(probs.mmdae; save_everystep = false)
        solve(probs.dae; save_everystep = false)
    end
end

end
