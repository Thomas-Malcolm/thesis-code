include("../../../src/ACSimulation.jl")


# Data 
current_df = DataFrame(CSV.File("../isttok_j_profile.csv"))
pressure_df = DataFrame(CSV.File("../isttok_p_profile.csv"))

current_data = Matrix(current_df)
pressure_data = Matrix(pressure_df)

# System Setup
cfg = Config(
    x0 = 0.46,      # m
    ρ = 0.085,      # m
)

val_h = 200

# Variations
a1_vals = range(-5.0, 5.0, length = val_h)
a2_vals = range(-20.0, 20.0, length = val_h)
a3_vals = range(-10.0, 10.0, length = val_h)

a1_expected = -0.04531
a2_expected = -1.0808
a3_expected = 2.1653


# Varying a1 and a2, with a fixed α
println("a1, a2, fixed α")

heat_vals = [
    sum_residuals(a1v, a2v, a3_expected)
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
savefig(hm1, "heat-a1-a2-large.pdf")

# Varying a1 and α, with a fixed a2
println("a1, a2, fixed α")

heat_vals = [
    sum_residuals(a1v, a2_expected, αv)
    for αv in x3_vals, a1v in x1_vals
]

hm2 = plot([
        heatmap(
            x = x1_vals,
            y = x3_vals,
            z = heat_vals,
        ),
    ],
    Layout(
        xaxis_title = "a1",
        yaxis_title = "alpha",
        showlegend = false
    )
)

savefig(hm2, "heat-a1-a3-large.pdf")

# Varying a2 and α, with a fixed a1
println("a2, α, fixed a1")
heat_vals = [
    sum_residuals(a1_expected, a2v, αv)
    for αv in x3_vals, a2v in x2_vals
]

hm3 = plot([
        heatmap(
            x = x2_vals,
            y = x3_vals,
            z = heat_vals,
        ),
        scatter(x = df1.x, y = df1.y),
        scatter(x = df2.x, y = df2.y)
    ],
    Layout(
        xaxis_title = "a2",
        yaxis_title = "alpha",
        showlegend = false
    )
)

savefig(hm3, "heat-a2-a3-large.pdf")