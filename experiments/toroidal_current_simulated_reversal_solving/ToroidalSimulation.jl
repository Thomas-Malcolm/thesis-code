include("../../src/ACSimulation.jl")

using Plots

using JuMP
using NLopt
using Ipopt


# Default
cfg = Config(
    B0 = 0.47
)
x_range = cfg.x_range
z_range = cfg.z_range
ρ = cfg.ρ


function current_inversion_system(initial_current_data::Data, initial_params::Parameters, h::Int = 50)
    a = Animation()

    M = maximum_current_for_params(cfg, initial_params)

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
        p = plot_for_start_params(previous_params, M)
        title!("$(i)/$(h)")

        frame(a, p)
    end

    gif(a, "simulated-inversion-fig-2-far.gif"; fps = 5)
end

function solve_for_parameters(current_data::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1v, a2v, αv)
        v = toroidal_residuals(cfg, Parameters(a1v, a2v, αv), current_data)
        # println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
        v
    end

    function ∇f(G::AbstractVector{T}, a1v::T, a2v::T, αv::T) where T
        l_params = Parameters(a1v, a2v, αv)

        g1,g2,g3 = residualToroidalDa1(cfg, l_params, current_data), residualToroidalDa2(cfg, l_params, current_data), residualToroidalDalpha(cfg, l_params, current_data)

        G[1] = g1
        G[2] = g2
        G[3] = g3
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
    # println("Optimising...")
    JuMP.optimize!(model)

    # Results
    println("\tObjective value: $(objective_value(model))")
    println("\ta1: $(value(a1))")
    println("\ta2: $(value(a2))")
    println("\tα : $(value(α))")
    # println("\tstopping reason: $(termination_status(model))")
    # println("\t..             : $(raw_status(model))")
    
    value(a1), value(a2), value(α)
end

function plot_for_start_params(params::Parameters, M::Float64)
    mag_field = magnetic_field_for_params(cfg, params)
    # Normalise current relative to the starting position
    curr_vals = current_profile_for_params_unnormalised(cfg, params) / M
    pres_vals = pressure_profile_for_params_normalised(cfg, params)

    p = contour(x_range, z_range, mag_field, xlims = (cfg.x0 - cfg.ρ, cfg.x0 + cfg.ρ), ylims = (-cfg.ρ, cfg.ρ))
    plot!(x_range, curr_vals, label = "Current")
    plot!(x_range, pres_vals, label = "Pressure")

    p
end



# Figure 1
# i_params = Parameters(-0.04531, -1.0808, 2.1683)
# ds = Data([(5.688442211055277, 0.23682196884837436), (5.447236180904523, 0.30576276376297173), (4.0201005025125625, -0.07973414793358713), (5.527638190954774, 0.28748266255678767), (4.71356783919598, 0.23248389606965747), (5.276381909547739, 0.3270127597649249), (5.567839195979899, 0.2765369764260688), (4.733668341708543, 0.23960739421881963), (5.396984924623116, 0.31458869633788256), (4.080402010050252, -0.052125947549691715), (5.457286432160804, 0.3037514301893421), (5.35678391959799, 0.3201259469155341), (4.341708542713568, 0.07204369953756981), (5.0954773869346734, 0.3225890543271238), (5.839195979899498, 0.17286041330620944), (5.417085427135678, 0.311307878206747), (5.236180904522613, 0.32836145903103364), (4.160804020100502, -0.014822850283991724), (5.768844221105527, 0.20463237273314155), (4.311557788944723, 0.05748246958014075), (4.954773869346734, 0.30149203683307585), (4.994974874371859, 0.3090397321598747), (5.507537688442211, 0.29251314225588265), (4.432160804020101, 0.1149692926137231), (5.045226130653266, 0.3167697487724618), (4.060301507537688, -0.06133800184445084), (5.075376884422111, 0.32049212442381647), (4.623115577889447, 0.19803804215363568), (4.522613065326633, 0.1558558406410243), (4.100502512562814, -0.0428806277570339), (5.256281407035176, 0.327859765896931), (5.1457286432160805, 0.32646630666345067), (5.989949748743719, 0.09693371451247428), (5.346733668341709, 0.32129361327193817), (5.698492462311558, 0.23304844923752568), (5.1959798994974875, 0.3283458584852187), (5.608040201005025, 0.2644379724026341), (4.422110552763819, 0.1102826609585354), (4.241206030150754, 0.023468126297322756), (4.231155778894473, 0.018633854354363516), (4.864321608040201, 0.28016326379411294), (5.577889447236181, 0.2736194549769976), (5.63819095477387, 0.25461620519055306), (4.78391959798995, 0.2564478005583263), (5.71859296482412, 0.22528616594625292), (4.8040201005025125, 0.2627673795599191), (4.301507537688442, 0.052618165622939125), ])

# Figure 2
i_params = Parameters(-0.01, 1.1, 4.566)
ds = Data([(5.939698492462312, 179338.5806059259), (4.562814070351759, -526936.5540833966), (4.442211055276382, -1.0854597284042102e6), (5.2462311557788945, 1.7872059255123525e6), (5.296482412060302, 1.472720574128773e6), (5.477386934673367, 305370.6877830321), (6.0, 363403.7865286629), (4.763819095477387, 1.0941245784127577e6), (4.592964824120603, -316437.2311204462), (4.1206030150753765, -195170.33400870944), (5.668341708542713, -281909.71702433855), (4.733668341708543, 842653.8712120074), (5.396984924623116, 794945.6112555265), (4.773869346733668, 1.1757970415618285e6), (4.080402010050252, 56721.25776864779), (5.115577889447236, 2.346659717310859e6), (5.467336683417085, 360787.49895301845), (4.924623115577889, 2.1498340798301157e6), (4.391959798994975, -1.1582990981108865e6), (4.341708542713568, -1.1349099155630243e6), (4.0703517587939695, 121311.33629234067), (5.768844221105527, -233289.85378274694), (4.492462311557789, -915010.4857137202), (4.311557788944723, -1.078062946001831e6), (4.14070351758794, -316167.0898646419), (5.035175879396985, 2.4232659956149794e6), (5.658291457286432, -274485.27797375136), (4.060301507537688, 186416.96366406046), (4.28140703517588, -992497.2226266149), (5.075376884422111, 2.41366058048281e6), (4.57286432160804, -459296.64836381347), (4.582914572864322, -389071.2545155596), (5.015075376884422, 2.4060521063843463e6), (4.482412060301508, -956699.4873253944), (4.522613065326633, -768396.2153420225), (5.21608040201005, 1.9552432452273981e6), (4.381909547738694, -1.1611362353221779e6), (5.376884422110553, 929286.8245819769), (4.100502512562814, -70659.74518812305), (4.894472361809045, 2.0031772541432146e6), (5.778894472361809, -217942.51264687825), (5.698492462311558, -289725.92291726393), (4.452261306532663, -1.059079286689106e6), (5.1959798994974875, 2.0560157878650536e6), (5.608040201005025, -198789.64087348946), (4.241206030150754, -839601.7654081613), (5.386934673366834, 861766.1116703277), (5.5376884422110555, 19351.571609320596), (5.105527638190955, 2.3686419731520903e6), ])

current_inversion_system(ds, i_params)