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

    gif(a, "simulated-inversion-fig-1.gif"; fps = 5)
end

function solve_for_parameters(current_data::Data, guess::Parameters)

    # Optimisation function and gradient
    function f(a1v, a2v, αv)
        v = un_residuals(cfg, Parameters(a1v, a2v, αv), current_data)
        # println("\tresiduals: $(v) $(round(a1v; digits = 5)), $(round(a2v; digits = 5)), $(round(αv; digits = 5))")
        v
    end

    function ∇f(G::AbstractVector{T}, a1v::T, a2v::T, αv::T) where T
        l_params = Parameters(a1v, a2v, αv)

        g1,g2,g3 = residualUnDa1(cfg, l_params, current_data), residualUnDa2(cfg, l_params, current_data), residualUnDalpha(cfg, l_params, current_data)

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
    # println("\tObjective value: $(objective_value(model))")
    # println("\ta1: $(value(a1))")
    # println("\ta2: $(value(a2))")
    # println("\tα : $(value(α))")
    # println("\tstopping reason: $(termination_status(model))")
    # println("\t..             : $(raw_status(model))")
    
    value(a1), value(a2), value(α)
end

function plot_for_start_params(params::Parameters, M::Float64)
    mag_field = magnetic_field_for_params(cfg, params)
    # Normalise current relative to the starting position
    curr_vals = current_profile_for_params_unnormalised(cfg, params) / M
    pres_vals = pressure_profile_for_params_normalised(cfg, params)

    p = contour(x_range, z_range, mag_field, xlims = (cfg.x0 - 1.0, cfg.x0 + 1.0), ylims = (-1.0, 1.0))
    plot!(x_range, curr_vals, label = "Current")
    plot!(x_range, pres_vals, label = "Pressure")

    p
end



# Figure 1
i_params = Parameters(-0.04531, -1.0808, 2.1683)
ds = Data([(5.688442211055277, 0.23682196884837436), (5.447236180904523, 0.30576276376297173), (4.0201005025125625, -0.07973414793358713), (5.527638190954774, 0.28748266255678767), (4.71356783919598, 0.23248389606965747), (5.276381909547739, 0.3270127597649249), (5.567839195979899, 0.2765369764260688), (4.733668341708543, 0.23960739421881963), (5.396984924623116, 0.31458869633788256), (4.080402010050252, -0.052125947549691715), (5.457286432160804, 0.3037514301893421), (5.35678391959799, 0.3201259469155341), (4.341708542713568, 0.07204369953756981), (5.0954773869346734, 0.3225890543271238), (5.839195979899498, 0.17286041330620944), (5.417085427135678, 0.311307878206747), (5.236180904522613, 0.32836145903103364), (4.160804020100502, -0.014822850283991724), (5.768844221105527, 0.20463237273314155), (4.311557788944723, 0.05748246958014075), (4.954773869346734, 0.30149203683307585), (4.994974874371859, 0.3090397321598747), (5.507537688442211, 0.29251314225588265), (4.432160804020101, 0.1149692926137231), (5.045226130653266, 0.3167697487724618), (4.060301507537688, -0.06133800184445084), (5.075376884422111, 0.32049212442381647), (4.623115577889447, 0.19803804215363568), (4.522613065326633, 0.1558558406410243), (4.100502512562814, -0.0428806277570339), (5.256281407035176, 0.327859765896931), (5.1457286432160805, 0.32646630666345067), (5.989949748743719, 0.09693371451247428), (5.346733668341709, 0.32129361327193817), (5.698492462311558, 0.23304844923752568), (5.1959798994974875, 0.3283458584852187), (5.608040201005025, 0.2644379724026341), (4.422110552763819, 0.1102826609585354), (4.241206030150754, 0.023468126297322756), (4.231155778894473, 0.018633854354363516), (4.864321608040201, 0.28016326379411294), (5.577889447236181, 0.2736194549769976), (5.63819095477387, 0.25461620519055306), (4.78391959798995, 0.2564478005583263), (5.71859296482412, 0.22528616594625292), (4.8040201005025125, 0.2627673795599191), (4.301507537688442, 0.052618165622939125), ])

# Figure 2
# i_params = Parameters(0.01, 3.1, 5.566)
# ds = Data([(5.969849246231155, 0.3403533223820269), (4.562814070351759, -0.662168003246996), (4.261306532663316, -1.1575511295373828), (5.2462311557788945, 2.245869203639298), (5.688442211055277, -0.36375057975434644), (5.527638190954774, 0.07670058952474423), (4.693467336683417, 0.6276819480711135), (5.798994974874372, -0.22957610397369366), (6.0, 0.45666666666666256), (4.763819095477387, 1.3749174958098913), (5.668341708542713, -0.3542581985845453), (4.743718592964824, 1.1653623542683582), (5.899497487437186, 0.07758371988923224), (4.080402010050252, 0.07127803472214478), (5.457286432160804, 0.52543966086227), (5.407035175879397, 0.9160882551991394), (4.663316582914573, 0.30575580077424314), (4.0703517587939695, 0.1524443212402583), (5.336683417085427, 1.5114262941517946), (5.788944723618091, -0.25264139897014176), (4.211055276381909, -0.8803188073094237), (4.64321608040201, 0.09578600492539502), (4.954773869346734, 2.84870682347035), (4.0, 0.7350000000000001), (5.507537688442211, 0.1907939690585076), (4.964824120603015, 2.889101193949353), (4.432160804020101, -1.3922393633645336), (4.060301507537688, 0.23425846554813573), (5.075376884422111, 3.0330953408127717), (4.381909547738694, -1.4591268274763385), (4.190954773869347, -0.7517540139360478), (4.522613065326633, -0.9655951625915259), (5.376884422110553, 1.1677762651095192), (5.175879396984925, 2.697342086361143), (5.346733668341709, 1.425338773619112), (4.532663316582915, -0.8956023002125234), (4.010050251256281, 0.6506725099206537), (4.180904522613066, -0.6843404721194117), (4.291457286432161, -1.2868082690139686), (4.914572864321608, 2.6441284884774245), (4.683417085427136, 0.5197957297519358), (5.577889447236181, -0.15251967001717515), (4.723618090452261, 0.9516646933735053), (4.824120603015075, 1.9604749501865035), (4.301507537688442, -1.3227026413117784), (4.844221105527638, 2.13578746422183), ])

current_inversion_system(ds, i_params)