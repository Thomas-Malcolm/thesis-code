include("../../src/ACSimulation.jl")

using Plots

using JuMP
using NLopt
using Ipopt


# Default
cfg = Config()
x_range = cfg.x_range
z_range = cfg.z_range


function current_inversion_system(initial_current_data::Data, initial_params::Parameters, h::Int = 50)
    a = Animation()

    # End state
    target_current = -initial_current_data

    # Track previous state's parameters to help next state's guess
    previous_params = initial_params

    # Iterate through current steps
    for i in 1:h
        println("[!] Current Step $(i)/$(h)")

        next_current::Data = Data()

        # For each data point, calculate its value at ith step
        for (j::Int, curr::Datum) in pairs(initial_current_data.d)
            sign = curr.val / abs(curr.val)

            new_d = Datum(
                curr.x, 
                curr.val - sign * (abs(target_current.d[j].val - curr.val) / h) * i
            )

            push!(next_current.d, new_d)
        end

        # Solve for parameter values given current data
        a1c, a2c, αc = solve_for_parameters(next_current, previous_params)
        previous_params = Parameters(a1c, a2c, αc)

        # Plot
        p = plot_for_start_params(previous_params)

        frame(a, p)
    end

    gif(a, "current-inversion-fig-1.gif"; fps = 5)
end

function solve_for_parameters(current_data::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1, a2, α)
        residuals(cfg, Parameters(a1, a2, α), current_data)
    end

    function ∇f(G::AbstractVector{T}, a1::T, a2::T, α::T) where T
        l_params = Parameters(a1, a2, α)

        G[1] = residualDa1(cfg, l_params, current_data)
        G[2] = residualDa2(cfg, l_params, current_data)
        G[3] = residualDalpha(cfg, l_params, current_data)
    end

    # Optimiser initialisation
    model = Model(NLopt.Optimizer)

    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    set_optimizer_attribute(model, "ftol_abs", 1e-5)

    register(model, :param_solve, 3, f, ∇f)

    # Tight search space around the previous values
    @variable(model, guess.a1 - 0.008 ≤ a1 ≤ guess.a1 + 0.008)
    @variable(model, guess.a2 - 0.2 ≤ a2 ≤ guess.a2 + 0.2)
    @variable(model, guess.α - 0.03 ≤ α ≤ guess.α + 0.03)

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
    println("\tα  : $(value(α))")
    
    value(a1), value(a2), value(α)
end

function plot_for_start_params(params::Parameters)
    mag_field = magnetic_field_for_params(cfg, params)
    curr_vals = current_profile_for_params_normalised(cfg, params)
    pres_vals = pressure_profile_for_params_normalised(cfg, params)

    p = contour(x_range, z_range, mag_field)
    plot!(x_range, curr_vals, label = "Current")
    plot!(x_range, pres_vals, label = "Pressure")

    p
end



# Figure 1
i_params = Parameters(-0.04531, -1.0808, 2.1683)
ds = Data([(4.562814070351759, 0.5271319912938993), (4.1708542713567835, -0.03071488391231543), (4.442211055276382, 0.3641464943325205), (4.874371859296482, 0.8608953285081333), (4.050251256281407, -0.200707612808032), (5.849246231155779, 0.5116083973086332), (4.693467336683417, 0.6853595225343578), (6.0, 0.2792101836871248), (4.633165829145729, 0.6150058837311638), (5.276381909547739, 0.9954062001816789), (4.733668341708543, 0.729349784351659), (4.743718592964824, 0.7399473971813753), (5.899497487437186, 0.43644828097097393), (4.080402010050252, -0.15866809423157385), (5.35678391959799, 0.9744431765531748), (5.467336683417085, 0.9182343507501194), (4.391959798994975, 0.29246728959561613), (4.0703517587939695, -0.17269812613700056), (4.311557788944723, 0.17497300919682354), (4.14070351758794, -0.07377749046144401), (4.954773869346734, 0.9177227304059346), (4.834170854271357, 0.8272421593617327), (5.557788944723618, 0.8504227415781456), (5.035175879396985, 0.9599844963148867), (5.165829145728643, 0.9967729184103393), (5.738693467336684, 0.6612569751687866), (5.809045226130653, 0.5687994533183357), (4.703517587939698, 0.6965887940507236), (5.819095477386934, 0.5547786010884644), (5.21608040201005, 1.0), (4.582914572864322, 0.5528281724139369), (5.517587939698492, 0.8828491946944586), (4.321608040201005, 0.18976790719521844), (5.376884422110553, 0.9665403309747135), (5.256281407035176, 0.9979844333857804), (5.1457286432160805, 0.9937428314320712), (5.778894472361809, 0.6096737484118069), (5.1959798994974875, 0.9994640685431111), (4.4623115577889445, 0.3922702205715146), (4.090452261306533, -0.14461259382704408), (4.231155778894473, 0.056720276514437375), (5.206030150753769, 0.9998591389659488), (5.869346733668341, 0.48197986791545677), (4.221105527638191, 0.04203981316528289)])


current_inversion_system(ds, i_params)