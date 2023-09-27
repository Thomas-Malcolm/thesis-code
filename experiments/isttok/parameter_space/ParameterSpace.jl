include("../../../src/ACSimulation.jl")

using PlotlyJS
using DataFrames
using CSV

# Data 
current_df = DataFrame(CSV.File("../isttok_j_profile.csv"))

curr_data = Matrix(current_df)[1,:] # Get first time slice




# System Setup
cfg = Config(
    x0 = 0.46,      # m
    ρ = 0.085,      # m
    B0 = 0.47       # T (Teslas)
)

val_h = 1000 # precision of mesh for heatmap

# Variations
a1_vals = range(-50.0, 50.0, length = val_h)
a2_vals = range(-50.0, 50.0, length = val_h)
a3_vals = range(-50.0, 50.0, length = val_h)

a1_expected = -0.05
a2_expected = -1.1
a3_expected = 2.17

# Data Collation

radius_vals = parse.(Float64, names(current_df)[2:end]) / 100 .+ cfg.x0 # Adjust to be centered, and in metres
curr_data = curr_data[2:end]
curr_d = Data(radius_vals, curr_data)


# Varying a1 and a2, with a fixed α
println("a1, a2, fixed α")

heat_vals = [
    toroidal_residuals(cfg, Parameters(a1v, a2v, a3_expected), curr_d)
    for a2v in a2_vals, a1v in a1_vals
]

hm1 = plot([    
        heatmap(
            x = a1_vals,
            y = a2_vals,
            z = heat_vals,
        ),
    ],
    Layout(
        xaxis_title = "a1",
        yaxis_title = "a2",
        showlegend = false
    )
)
savefig(hm1, "heat-a1-a2-gigantic.pdf")

# Varying a1 and α, with a fixed a2
println("a1, a3, fixed α")

heat_vals = [
    toroidal_residuals(cfg, Parameters(a1v, a2_expected, αv), curr_d)
    for αv in a3_vals, a1v in a1_vals
]

hm2 = plot([
        heatmap(
            x = a1_vals,
            y = a3_vals,
            z = heat_vals,
        ),
    ],
    Layout(
        xaxis_title = "a1",
        yaxis_title = "alpha",
        showlegend = false
    )
)

savefig(hm2, "heat-a1-a3-gigantic.pdf")

# Varying a2 and α, with a fixed a1
println("a2, α, fixed a1")
heat_vals = [
    toroidal_residuals(cfg, Parameters(a1_expected, a2v, αv), curr_d)
    for αv in a3_vals, a2v in a2_vals
]

hm3 = plot([
        heatmap(
            x = a2_vals,
            y = a3_vals,
            z = heat_vals,
        ),
    ],
    Layout(
        xaxis_title = "a2",
        yaxis_title = "alpha",
        showlegend = false
    )
)

savefig(hm3, "heat-a2-a3-gigantic.pdf")