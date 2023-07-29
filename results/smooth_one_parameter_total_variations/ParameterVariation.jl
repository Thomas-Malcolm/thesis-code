include("../../src/ACSimulation.jl")
using .ACSimulation
using Plots


# Default
cfg = Config()
x_range = cfg.x_range
z_range = cfg.z_range

exp_params = Parameters(0.01, 3.1, 5.566)
frame_count = 200

println("Creating animation a1")
anim = Animation()
i = 0

for a1 in range(-0.06, 0.06, length = frame_count)
    global i
    println("$(i) / $(frame_count)")

    params = Parameters(a1, exp_params.a2, exp_params.α)

    mag_field = magnetic_field_for_params(cfg, params)
    curr_vals = current_profile_for_params_normalised(cfg, params)
    pres_vals = pressure_profile_for_params_normalised(cfg, params)


    if (abs(a1 - 0.01) ≤ 1e-5 )
        lcol = :red
    else
        lcol = :black
    end

    contour(x_range, z_range, mag_field, title = "a1 = $(round(a1; digits = 3)) ($(i) / $(frame_count))")
    plot!(x_range, curr_vals, linecolor = lcol)
    plot!(x_range, pres_vals, linecolor = lcol)

    i += 1

    frame(anim)
end

gif(anim, "total-a1-fig-1.gif" ; fps = 30)



println("Creating animation a2")
anim = Animation()
i = 0

for a2 in range(2.0, 4.0, length = frame_count)
    global i
    println("$(i) / $(frame_count)")

    params = Parameters(exp_params.a1, a2, exp_params.α)

    mag_field = magnetic_field_for_params(cfg, params)
    curr_vals = current_profile_for_params_normalised(cfg, params)
    pres_vals = pressure_profile_for_params_normalised(cfg, params)

    if (abs(a2 - 3.1) ≤ 1e-4 )
        lcol = :red
    else
        lcol = :blue
    end

    contour(x_range, z_range, mag_field, title = "a2 = $(round(a2; digits = 3)) ($(i) / $(frame_count))")
    plot!(x_range, curr_vals, linecolor = lcol)
    plot!(x_range, pres_vals, linecolor = :purple)

    i += 1

    frame(anim)
end

gif(anim, "total-a2-fig-1.gif" ; fps = 30)



println("Creating animation α")
anim = Animation()
i = 0

for α in range(4.8, 6.8, length = frame_count)
    global i
    println("$(i) / $(frame_count)")

    mag_field = magnetic_field_for_params(cfg, params)
    curr_vals = current_profile_for_params_normalised(cfg, params)
    pres_vals = pressure_profile_for_params_normalised(cfg, params)

    if (abs(α - 5.566) ≤ 1e-4 )
        lcol = :red
    else
        lcol = :black
    end

    contour(x_range, z_range, mag_field, ylims = (-5, 5), title = "a2 = $(round(a2; digits = 3)) ($(i) / $(frame_count))")
    plot!(x_range, curr_vals, linecolor = lcol)
    plot!(x_range, pres_vals, linecolor = :purple)

    i += 1

    frame(anim)
end

gif(anim, "total-alpha-fig-1.gif" ; fps = 10)