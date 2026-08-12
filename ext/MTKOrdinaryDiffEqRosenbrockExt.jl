module MTKOrdinaryDiffEqRosenbrockExt

using ModelingToolkit
using ModelingToolkitBase
using OrdinaryDiffEqRosenbrock: Rodas5P
using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    probs = ModelingToolkitBase._solver_precompile_problems()
    @compile_workload begin
        ModelingToolkitBase._solve_precompile(probs.ode, Rodas5P())
        ModelingToolkitBase._solve_precompile(probs.mmdae, Rodas5P())
    end
end

end
