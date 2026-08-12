module MTKOrdinaryDiffEqTsit5Ext

using DiffEqBase
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits, D_nounits
using OrdinaryDiffEqTsit5: Tsit5
using PrecompileTools: @compile_workload, @setup_workload

# Tsit5 precompiles handwritten AutoDespecialize (`sys = nothing`). MTK problems
# keep `sys::System` and an ObservedFunctionCache in the ODEFunction type.
@setup_workload begin
    @parameters a = 1.0 b = 1.0
    @variables x(t_nounits) y(t_nounits)
    prob = ODEProblem(
        mtkcompile(
            System(
                [D_nounits(x) ~ a * y, D_nounits(y) ~ -b * x],
                t_nounits; name = :precompile_tsit5
            )
        ),
        [x => 1.0, y => 0.0], (0.0, 1.0)
    )
    @compile_workload begin
        solve(prob, Tsit5(); adaptive = false, dt = 0.1, save_everystep = false)
        wrapped = DiffEqBase.get_concrete_problem(prob, true; alg = Tsit5())
        solve(wrapped, Tsit5(); adaptive = false, dt = 0.1, save_everystep = false)
    end
end

end
