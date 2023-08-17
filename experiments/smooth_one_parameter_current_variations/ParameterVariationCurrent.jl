include("../../src/ACSimulation.jl")
using Plots


# Default
cfg = Config()
x_range = cfg.x_range

exp_params = Parameters(0.01, 3.1, 5.566)

frame_count = 200


println("Creating animation a1")
anim = Animation()
i = 0

for a1 in range(-0.06, 0.06, length = frame_count)
    global i
    println("$(i) / $(frame_count)")

    params = Parameters(a1, exp_params.a2, exp_params.α)

    curr_vals = current_profile_for_params_normalised(cfg, params)

    if (abs(a1 - 0.01) ≤ 1e-4 )
        lcol = :red
    else
        lcol = :black
    end

    plot(x_range, curr_vals, linecolor = lcol, title = "a1 = $(a1) ($(i) / $(frame_count))")

    i += 1

    frame(anim)
end

gif(anim, "current-a1-variation.gif" ; fps = 30)


println("Creating animation a2")
anim = Animation()
i = 0

for a2 in range(2.0, 4.0, length = frame_count)
    global i
    println("$(i) / $(frame_count)")

    params = Parameters(exp_params.a1, a2, exp_params.α)

    curr_vals = current_profile_for_params_normalised(cfg, params)

    if (abs(a2 - 3.1) ≤ 1e-4 )
        lcol = :red
    else
        lcol = :black
    end

    plot(x_range, curr_vals, linecolor = lcol, title = "a2 = $(a2) ($(i) / $(frame_count))")

    i += 1

    frame(anim)
end

gif(anim, "current-a2-variation.gif" ; fps = 30)



println("Creating animation α")
anim = Animation()
i = 0

for α in range(4.8, 6.8, length = frame_count)
    global i
    println("$(i) / $(frame_count)")

    params = Parameters(exp_params.a1, exp_params.a2, α)

    curr_vals = current_profile_for_params_normalised(cfg, params)

    if (abs(α - 5.566) ≤ 1e-4 )
        lcol = :red
    else
        lcol = :black
    end

    plot(x_range, curr_vals, linecolor = lcol, ylims = (-10, 10), title = "α = $(α) ($(i) / $(frame_count))")

    i += 1

    frame(anim)
end

gif(anim, "current-alpha-variation.gif" ; fps = 10)