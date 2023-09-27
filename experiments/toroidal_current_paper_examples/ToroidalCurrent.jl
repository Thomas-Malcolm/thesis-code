include("../../src/ACSimulation.jl")

using Plots


 # Setup
cfg = Config(
    x0 = 5.0,      # m
    ρ = 1.0,      # m
    B0 = 0.47       # T
)

params = Parameters(
    0.01, 3.1, 5.566
)


x_range = cfg.x_range
z_range = cfg.z_range


# Toroidal Current
tcurr = toroidal_current_density_profile(cfg, params)
p1 = plot(x_range, tcurr)
savefig(p1, "toroidal-current.pdf")

# Normalised current wrt reactor
mcurr = current_profile_for_params_unnormalised(cfg, params)
p2 = plot(x_range, mcurr)
savefig(p2, "modified-current.pdf")

println("Tcurr:")
println(tcurr)
println("Mcurr:")
println(mcurr)