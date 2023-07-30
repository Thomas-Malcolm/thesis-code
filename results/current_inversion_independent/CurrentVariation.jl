include("../../src/ACSimulation.jl")
using Plots

# Default
cfg = Config(h = 400)
x_range = cfg.x_range
z_range = cfg.z_range

# Figure 1
# params = Parameters(-0.04531, -1.0808, 2.1683)
# Figure 2
params = Parameters(0.01, 3.1, 5.566)


function smoothly_invert(curr::Vector{Float64}, h::Int = 200)
    a = Animation()
    M = maximum(abs.(curr))

    targetCurrent = -curr

    for i in 1:h
        println("$(i)/$(h)")

        newCurrent::Vector{Float64} = Vector{Float64}()

        for j in 1:length(curr)
            if curr[j] > 0
                push!(
                    newCurrent,
                    curr[j] - (abs(targetCurrent[j] - curr[j]) / h) * i
                )
            else
                push!(
                    newCurrent,
                    curr[j] + (abs(targetCurrent[j] - curr[j]) / h) * i
                )
            end
        end

        p = scatter(
            x_range, 
            newCurrent,
            title = "Step $(i)",
            xlabel = "x",
            ylabel = "Current Vals",
            xlims = (4.0, 6.0),
            ylims = (-M, M),
            msw = 0 # marker stroke width
        )

        frame(a)
    end

    gif(a, "current-inversion.gif"; fps = 30)
end

curr_vals = current_profile_for_params_normalised(cfg, params)

smoothly_invert(curr_vals)

