module MTKOrdinaryDiffEqTsit5Ext

using ModelingToolkitBase
using OrdinaryDiffEqTsit5: Tsit5
using PrecompileTools: @compile_workload, @setup_workload

# Tsit5 precompiles handwritten AutoDespecialize (`sys = nothing`). MTK problems
# keep `sys::System` and an ObservedFunctionCache in the ODEFunction type.
@setup_workload begin
    probs = ModelingToolkitBase._solver_precompile_problems()
    @compile_workload begin
        ModelingToolkitBase._solve_precompile(probs.ode, Tsit5())
    end
end

end
