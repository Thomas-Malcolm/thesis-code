include("../../../src/ACSimulation.jl")

using Plots

using CSV
using DataFrames

using JuMP
using NLopt

# Data 
current_df = DataFrame(CSV.File("../isttok_j_profile.csv"))
pressure_df = DataFrame(CSV.File("../isttok_p_profile.csv"))

current_data = Matrix(current_df)
pressure_data = Matrix(pressure_df)

# System Setup
cfg = Config(
    x0 = 46.0,      # cm
    ρ = 8.5,        # cm
)

x_range = cfg.x_range
z_range = cfg.z_range

radius_vals = parse.(Float64, names(current_df)[2:end]) .+ cfg.x0 # Adjust to be centered



# Helper Functions
function solve_for_parameters(cdata::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1v, a2v, αv)
        v = un_residuals(cfg, Parameters(a1v, a2v, αv), cdata)
        println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
        v
    end

    function ∇f(G::AbstractVector{T}, a1v::T, a2v::T, αv::T) where T
        l_params = Parameters(a1v, a2v, αv)

        G[1] = residualUnDa1(cfg, l_params, cdata)
        G[2] = residualUnDa2(cfg, l_params, cdata)
        G[3] = residualUnDalpha(cfg, l_params, cdata)
    end

    # Optimiser initialisation
    model = Model(NLopt.Optimizer)

    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    set_optimizer_attribute(model, "ftol_abs", 1e-7)
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
    JuMP.optimize!(model)

    # Optimisation Results
    println("[!] Optimisation Results")
    println("\tObjective value: $(objective_value(model))")
    println("\ta1: $(value(a1))")
    println("\ta2: $(value(a2))")
    println("\tα : $(value(α))")
    println("\tstopping reason: $(termination_status(model))")

    Parameters(value(a1), value(a2), value(α))
end





# "main"
i = 0
guess_params = Parameters(0.1, 0.1, 0.1) # Random guess to start off
for (curr_data, pres_data) in zip(eachrow(current_data), eachrow(pressure_data))
    global i
    global guess_params
    println("[+] Current $(i) / $(size(curr_data, 1))")

    # Data Collation
    time_slice = pres_data[1]
    
    curr_data = curr_data[2:end]
    pres_data = pres_data[2:end]

    curr_data_norm = curr_data / maximum(curr_data)
    pres_data_norm = pres_data / maximum(pres_data)

    ## `Data` structure 
    curr_d = Data(radius_vals, curr_data)
    pres_d = Data(radius_vals, pres_data)

    # Solve for `a1`, `a2`, and `alpha` here
    solved_params = solve_for_parameters(curr_d, guess_params)

    # Plotting
    mag_field = magnetic_field_for_params(cfg, solved_params)
    p = contour(x_range, z_range, mag_field, xlims = (cfg.x0 - cfg.ρ, cfg.x0 + cfg.ρ), ylims = (-cfg.k, cfg.k))
    plot!(radius_vals, curr_data_norm .* cfg.k, label = "Current")
    plot!(radius_vals, pres_data_norm .* cfg.k, label = "Pressure")
    
    title!("Time $(time_slice)")
    xlabel!("x (cm)")
    savefig(p, "graphs/slice_$(i).pdf")

    ## Comparison with simulated current / pressure
    p1 = plot(radius_vals, curr_data_norm, label = "Ip (Data)")
    plot!(radius_vals, pres_data_norm, label = "P (Data)")
    sim_curr = current_profile_for_params_normalised(cfg, solved_params) # / firstM
    sim_pres = pressure_profile_for_params_normalised(cfg, solved_params)
    p2 = plot(x_range, sim_curr, label = "Ip (Sim)")
    plot!(x_range, sim_pres, label = "P (Sim)")

    p = plot(
        p1, p2, layout = (2,1)
    )
    
    savefig(p, "graphs/comparison_$(i).pdf")

    # Continuation
    guess_params = solved_params
    i += 1
end