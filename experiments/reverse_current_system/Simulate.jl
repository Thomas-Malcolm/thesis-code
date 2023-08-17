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
        title!("$(i)/$(h)")

        frame(a, p)
    end

    gif(a, "current-inversion-fig-2.gif"; fps = 5)
end

function solve_for_parameters(current_data::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1v, a2v, αv)
        v = residuals(cfg, Parameters(a1v, a2v, αv), current_data)
        println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
        println("banana")
        v
    end

    function ∇f(G::AbstractVector{T}, a1v::T, a2v::T, αv::T) where T
        l_params = Parameters(a1v, a2v, αv)

        g1,g2,g3 = residualDa1(cfg, l_params, current_data), residualDa2(cfg, l_params, current_data), residualDalpha(cfg, l_params, current_data)

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
    println("Optimising...")
    JuMP.optimize!(model)

    # Results
    println("\tObjective value: $(objective_value(model))")
    println("\ta1: $(value(a1))")
    println("\ta2: $(value(a2))")
    println("\tα : $(value(α))")
    println("\tstopping reason: $(termination_status(model))")
    println("\t..             : $(raw_status(model))")
    
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
# i_params = Parameters(-0.04531, -1.0808, 2.1683)
# ds = Data([(4.562814070351759, 0.5271319912938993), (4.1708542713567835, -0.03071488391231543), (4.442211055276382, 0.3641464943325205), (4.874371859296482, 0.8608953285081333), (4.050251256281407, -0.200707612808032), (5.849246231155779, 0.5116083973086332), (4.693467336683417, 0.6853595225343578), (6.0, 0.2792101836871248), (4.633165829145729, 0.6150058837311638), (5.276381909547739, 0.9954062001816789), (4.733668341708543, 0.729349784351659), (4.743718592964824, 0.7399473971813753), (5.899497487437186, 0.43644828097097393), (4.080402010050252, -0.15866809423157385), (5.35678391959799, 0.9744431765531748), (5.467336683417085, 0.9182343507501194), (4.391959798994975, 0.29246728959561613), (4.0703517587939695, -0.17269812613700056), (4.311557788944723, 0.17497300919682354), (4.14070351758794, -0.07377749046144401), (4.954773869346734, 0.9177227304059346), (4.834170854271357, 0.8272421593617327), (5.557788944723618, 0.8504227415781456), (5.035175879396985, 0.9599844963148867), (5.165829145728643, 0.9967729184103393), (5.738693467336684, 0.6612569751687866), (5.809045226130653, 0.5687994533183357), (4.703517587939698, 0.6965887940507236), (5.819095477386934, 0.5547786010884644), (5.21608040201005, 1.0), (4.582914572864322, 0.5528281724139369), (5.517587939698492, 0.8828491946944586), (4.321608040201005, 0.18976790719521844), (5.376884422110553, 0.9665403309747135), (5.256281407035176, 0.9979844333857804), (5.1457286432160805, 0.9937428314320712), (5.778894472361809, 0.6096737484118069), (5.1959798994974875, 0.9994640685431111), (4.4623115577889445, 0.3922702205715146), (4.090452261306533, -0.14461259382704408), (4.231155778894473, 0.056720276514437375), (5.206030150753769, 0.9998591389659488), (5.869346733668341, 0.48197986791545677), (4.221105527638191, 0.04203981316528289)])

# Figure 2
i_params = Parameters(0.02, 3.2, 5.566)
ds = Data([(5.939698492462312, 0.07391283429660968), (4.261306532663316, -0.37964391815023446), (5.085427135678392, 0.9900437315215419), (4.874371859296482, 0.7789893924915741), (5.527638190954774, 0.025155616532674587), (4.71356783919598, 0.2767671253314692), (6.0, 0.1497737060575689), (4.763819095477387, 0.45093391723541876), (5.005025125628141, 0.9858129071187706), (4.613065326633166, -0.06787359657945002), (4.733668341708543, 0.3472924550972594), (5.065326633165829, 0.9980120975882734), (5.396984924623116, 0.3276299112049809), (4.542713567839196, -0.2694549538699903), (4.0703517587939695, 0.04999741086040349), (4.391959798994975, -0.47738288669021317), (5.236180904522613, 0.7605113200751149), (4.814070351758794, 0.6128513018709567), (5.437185929648241, 0.2217655419150166), (5.768844221105527, -0.09614838172280431), (4.64321608040201, 0.03141509112290166), (4.0, 0.24105914004156384), (5.306532663316583, 0.5795581749555674), (4.964824120603015, 0.9475431963353556), (4.432160804020101, -0.45661499818322276), (4.060301507537688, 0.07683012823468013), (4.28140703517588, -0.40904908753037333), (5.819095477386934, -0.05851199980513341), (5.185929648241206, 0.8665683299568029), (4.522613065326633, -0.31668780887429554), (4.150753768844221, -0.1546029804574642), (5.366834170854271, 0.4110405302559888), (4.894472361809045, 0.8255920614068505), (4.090452261306533, -0.003005069711790827), (4.4623115577889445, -0.4240091975918829), (5.8793969849246235, 0.0024705599975605465), (4.884422110552764, 0.8028949476662932), (5.577889447236181, -0.050022123120766426), (4.944723618090452, 0.9196113760569193), (4.221105527638191, -0.30867021059354893), (4.824120603015075, 0.6429801436258316), (4.8040201005025125, 0.5818830929634634), (5.889447236180905, 0.013829231772188146), (4.673366834170854, 0.13525518964728409), (5.105527638190955, 0.976215162839312), ])

current_inversion_system(ds, i_params)