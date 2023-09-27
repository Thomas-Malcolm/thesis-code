include("../../../src/ACSimulation.jl")

using Plots

using DataFrames
using CSV

using JuMP
using NLopt

# Data
current_df = DataFrame(CSV.File("../isttok_j_profile_slice_1_carved.csv"))
radius_vals = parse.(Float64, names(current_df)) / 100 

c_data = Float64.(Matrix(current_df)[1,:])

# System Setup
cfg = Config(
    x0 = 0.46,      # m
    ρ = 0.085,      # m
    B0 = 0.47       # T
)

x_range = cfg.x_range
z_range = cfg.z_range


# Scale data
radius_vals = parse.(Float64, names(current_df)) / 100 .+ cfg.x0 # Adjust to be centered, and in metres

B0 = cfg.B0
μ0 = cfg.μ0
a = cfg.ρ

scaled_current = c_data .* ((μ0 * a) / B0)

curr_d = Data(radius_vals, scaled_current)

function solve_for_parameters(current_data::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1v, a2v, αv)
        v = un_residuals(cfg, Parameters(a1v, a2v, αv), current_data)
        println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
        println("banana")
        v
    end

    function ∇f(G::AbstractVector{T}, a1v::T, a2v::T, αv::T) where T
        l_params = Parameters(a1v, a2v, αv)

        g1,g2,g3 = residualUnDa1(cfg, l_params, current_data), residualUnDa2(cfg, l_params, current_data), residualUnDalpha(cfg, l_params, current_data)

        println("\t\tdfa1: $(g1)")
        println("\t\tdfa2: $(g2)")
        println("\t\tdfal: $(g3)")

        G[1] = g1
        G[2] = g2
        G[3] = g3

        # G[1] = residualDa1(cfg, l_params, current_data)
        # G[2] = residualDa2(cfg, l_params, current_data)
        # G[3] = residualDalpha(cfg, l_params, current_data)
    end

    # Optimiser initialisation
    model = Model(NLopt.Optimizer)

    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    set_optimizer_attribute(model, "ftol_abs", 1e-10)
    set_optimizer_attribute(model, "xtol_rel", 1e-32)
    set_optimizer_attribute(model, "verbosity", 1)

    register(model, :param_solve, 3, f, ∇f)

    # Tight search space around the previous values
    # @variable(model, guess.a1 - 0.05 ≤ a1 ≤ guess.a1 + 0.05)
    # @variable(model, guess.a2 - 0.3 ≤ a2 ≤ guess.a2 + 0.3)
    # @variable(model, guess.α - 0.08 ≤ α ≤ guess.α + 0.08)
    @variable(model, a1)
    @variable(model, a2)
    @variable(model, α)

    set_start_value(a1, guess.a1)
    set_start_value(a2, guess.a2)
    set_start_value(α, guess.α)

    @NLobjective(model, Min, param_solve(a1, a2, α))

    # Run optimisation
    println("Optimising...")
    JuMP.optimize!(model)

    # Results
    println("\tObjective value: $(objective_value(model))")
    println("\ta1: $(value(a1))")
    println("\ta2: $(value(a2))")
    println("\tα : $(value(α))")
    println("\tstopping reason: $(termination_status(model))")
    println("\t..             : $(raw_status(model))")
    
    Parameters(value(a1), value(a2), value(α))
end

first_guess = Parameters(-5, 5.0, 20.0)

params = solve_for_parameters(curr_d, first_guess)
p1 = plot(radius_vals, scaled_current, label = "Ip (Data)")

# sim_curr = current_profile_for_params_unnormalised(cfg, params)
sim_curr = [ current_profile_for_params_unnormalised_x(xx, cfg, params) for xx in radius_vals ]
p2 = plot(radius_vals, sim_curr, label = "Ip (Sim)")

p = plot(
    p1, p2, layout = (2,1)
)

savefig(p, "graphs/comparison_reverse_scale_0.pdf") 