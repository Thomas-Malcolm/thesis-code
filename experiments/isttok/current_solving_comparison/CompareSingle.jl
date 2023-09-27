include("../../../src/ACSimulation.jl")

using Plots

using CSV
using DataFrames

using JuMP
using NLopt

# Data
current_df = DataFrame(CSV.File("../isttok_j_profile_slice_1_carved.csv"))
radius_vals = parse.(Float64, names(current_df)) / 100 

current_data = Matrix(current_df)

# System Setup
cfg = Config(
    x0 = 0.46,      # m
    ρ = 0.085,      # m
    B0 = 0.47       # T
)

x_range = cfg.x_range
z_range = cfg.z_range

radius_vals = parse.(Float64, names(current_df)) / 100 .+ cfg.x0 # Adjust to be centered, and in metres


fixed_alpha = 5.5

# Helper Functions
function solve_for_parameters(cdata::Data, guess::Parameters)

    # Optimisation function and gradient
    # function f(a1v, a2v, αv)
    #     v = toroidal_residuals(cfg, Parameters(a1v, a2v, αv), cdata)
    #     # println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
    #     v
    # end

    function f(a1v, a2v)
        v = toroidal_residuals(cfg, Parameters(a1v, a2v, 5.5), cdata)
        # println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
        v
    end

    # function ∇f(G::AbstractVector{T}, a1v::T, a2v::T, αv::T) where T
    #     l_params = Parameters(a1v, a2v, αv)

    #     g1 = residualToroidalDa1(cfg, l_params, cdata)
    #     println("Da1: $(g1)")
    #     g2 = residualToroidalDa2(cfg, l_params, cdata)
    #     println("Da2: $(g2)")
    #     g3 = residualToroidalDalpha(cfg, l_params, cdata)
    #     println("Da3: $(g3)")


    #     G[1] = g1
    #     G[2] = g2
    #     G[3] = g3
    # end

    function ∇f(G::AbstractVector{T}, a1v::T, a2v::T) where T
        l_params = Parameters(a1v, a2v, 5.5)

        g1 = residualToroidalDa1(cfg, l_params, cdata)
        println("Da1: $(g1)")
        g2 = residualToroidalDa2(cfg, l_params, cdata)
        println("Da2: $(g2)")


        G[1] = g1
        G[2] = g2
    end

    # Optimiser initialisation
    model = Model(NLopt.Optimizer)

    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    set_optimizer_attribute(model, "ftol_abs", 1e-7)
    set_optimizer_attribute(model, "xtol_rel", 1e-32)
    # set_optimizer_attribute(model, "verbosity", 1)

    register(model, :param_solve, 2, f, ∇f)

    # Tight search space around the previous values
    # @variable(model, guess.a1 - 0.05 ≤ a1 ≤ guess.a1 + 0.05)
    # @variable(model, guess.a2 - 0.3 ≤ a2 ≤ guess.a2 + 0.3)
    # @variable(model, guess.α - 0.08 ≤ α ≤ guess.α + 0.08)
    @variable(model, a1)
    @variable(model, a2)
    # @variable(model, α)

    set_start_value(a1, guess.a1)
    set_start_value(a2, guess.a2)
    # set_start_value(α, guess.α)

    @NLobjective(model, Min, param_solve(a1, a2))#, α))

    # Run optimisation
    JuMP.optimize!(model)

    # Optimisation Results
    println("[!] Optimisation Results")
    println("\tObjective value: $(objective_value(model))")
    println("\ta1: $(value(a1))")
    println("\ta2: $(value(a2))")
    # println("\tα : $(value(α))")
    println("\tstopping reason: $(termination_status(model))")

    Parameters(value(a1), value(a2), 5.5)
end

curr_data = Float64.(current_data[1,:])
curr_d = Data(radius_vals, curr_data)

guess_params = Parameters(-0.1, 2.0, 5.1)
solved_params = solve_for_parameters(curr_d, guess_params)

p1 = plot(radius_vals, curr_data, label = "Ip (Data)")
sim_curr = toroidal_current_density_profile(cfg, solved_params) 
p2 = plot(x_range, sim_curr, label = "Ip (Sim)")

p = plot(
    p1, p2, layout = (2,1)
)

savefig(p, "graphs/comparison_edited_a1_a2_0.pdf")