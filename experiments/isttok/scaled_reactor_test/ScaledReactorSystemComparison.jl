include("../../../src/ACSimulation.jl")

using Plots

using DataFrames
using CSV

using JuMP
using NLopt

# Data
current_df = DataFrame(CSV.File("../isttok_j_profile.csv"))
pressure_df = DataFrame(CSV.File("../isttok_p_profile.csv"))

current_data = Matrix(current_df)
pressure_data = Matrix(pressure_df)

# System Setup
cfg = Config(
    x0 = 5.0,       # 0.46,      # m
    ρ = 0.923913,   # 0.085,          # m
    B0 = 0.47           # T
)

origx0 = 0.46
origr = 0.085
x_range = cfg.x_range
z_range = cfg.z_range
x_plot_range = range(cfg.x0 - origr, cfg.x0 + origr, length = cfg.h)
z_plot_range = range(-origr, origr, length = cfg.h)

radius_vals = parse.(Float64, names(current_df)[2:end]) / 100 .+ cfg.x0 # Adjust to be centered, and in metres


# Helper Functions
function solve_for_parameters(current_data::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1v, a2v, αv)
        v = un_residuals(cfg, Parameters(a1v, a2v, αv), current_data)
        println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
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
    set_optimizer_attribute(model, "ftol_abs", 1e-18)
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




# "main"
i = 0
guess_params = Parameters(-5, 5.0, 20.0) # Random guess to start off
first_guess = Parameters(-5, 5.0, 20.0) # Use a single guess
for (curr_data, pres_data) in zip(eachrow(current_data), eachrow(pressure_data))
    global i
    global guess_params
    B0 = cfg.B0
    μ0 = cfg.μ0
    a = cfg.ρ

    println("[+] Current $(i) / $(size(curr_data, 1))")

    # Data Collation
    time_slice = pres_data[1]
    curr_data = Float64.(curr_data[2:end])
    pres_data = Float64.(pres_data[2:end])

    scaled_current = curr_data .* ((μ0 * a) / B0)
    scaled_pres = pres_data .* ((2 * μ0) / B0^2)

    ## `Data` structure 
    curr_d = Data(radius_vals, scaled_current)
    println(curr_d)
    
    # Solve for `a1`, `a2`, and `alpha` here
    solved_params = solve_for_parameters(curr_d, first_guess)

    ## Simulated system
    mag_field = magnetic_field_for_params_and_range(cfg, solved_params, x_plot_range, z_plot_range)
    sim_curr = [ current_profile_for_params_unnormalised_x(xx, cfg, solved_params) for xx in x_plot_range ]
    sim_pressure = [ pressure_profile_for_params_unnormalised_x(xx, cfg, solved_params) for xx in x_plot_range ]
    p1 = contour(x_plot_range, z_plot_range, mag_field, xlims = (cfg.x0 - origr, cfg.x0 + origr), ylims = (-origr, origr), legend = false)

    # Current
    p2 = plot(x_plot_range, sim_curr, label = "Sim", color = :blue, legend = :outerright, title = "Current Profile $(time_slice)")
    plot!(radius_vals, scaled_current, label = "Data", color = :black)

    # Pressure
    p3 = plot(x_plot_range, sim_pressure, color = :blue, legend = false, title = "Pressure Profile")
    plot!(radius_vals, scaled_pres, color = :black, legend = false)
    

    p = plot(
        p2, p3, p1, layout = @layout [a b ; c ]
    )
    
    savefig(p, "graphs/system_comparison_$(i).pdf")

    # Continuation
    guess_params = solved_params
    i += 1
end