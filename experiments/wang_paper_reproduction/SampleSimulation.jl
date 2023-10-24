include("../../src/ACSimulation.jl")
using Plots


# Default
cfg = Config()
x_range = cfg.x_range
z_range = cfg.z_range

# Figure 1
# params = Parameters(-0.04531, -1.0808, 2.1683)
# Figure 2
#params = Parameters(0.01, 3.1, 5.566)
# Figure 3
params = Parameters(-0.06, -0.0006, 5.50)

# Magnetic Field Lines
mag_field = magnetic_field_for_params(cfg, params)
contour(x_range, z_range, mag_field, title = "Wang Figure 3", xlabel = "x", ylabel = "z,j,p")

# Current Profile
curr_prof = current_profile_for_params_normalised(cfg, params)
plot!(x_range, curr_prof, label = "j")

# Pressure Profile
pres_prof = pressure_profile_for_params_normalised(cfg, params)
plot!(x_range, pres_prof, label = "p")


savefig("imgs/figure-3.pdf")
