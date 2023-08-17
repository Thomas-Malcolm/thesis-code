include("../../src/ACSimulation.jl")

"""
This will generate `data_size` randomly selected data points 
    for a given parameter set, and print a Julia style array
    of tuples of (x, value). 

The internal `Datum` structure can parse these and store them 
    for processing.
"""

cfg = Config()
x_range = cfg.x_range

# Fig 1
params = Parameters(-0.04531, -1.0808, 2.1683)

# Fig 2
params = Parameters(0.01, 3.1, 5.566)

# Generate Data
data_size = 50

x_vals = Set([rand(x_range) for _ in 1:data_size])

println("Current:")
print("[")
for x in x_vals
    print("($(x), $(current_profile_for_params_normalised_x(x, cfg, params))), ")
end
print("]")

## Need these for buffering
println()
exit(1)