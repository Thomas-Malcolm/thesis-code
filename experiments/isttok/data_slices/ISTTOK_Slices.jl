using CSV
using DataFrames

using Plots

# Current Data
# current_df = DataFrame(CSV.File("../isttok_j_profile.csv"))
# radius_vals = parse.(Float64, names(current_df)[2:end]) / 100 # convert to metres

# current_data = Matrix(current_df)


###########
current_df = DataFrame(CSV.File("../isttok_j_profile_slice_1_carved.csv"))
radius_vals = parse.(Float64, names(current_df)) / 100 

current_data = Matrix(current_df)


# Unnormalised
for curr_data in eachrow(current_data)
    p = plot(
        radius_vals, curr_data, 
        title = "ISTTOK Time Slice 1",
        xlabel = "radius (m)", ylabel = "Current Density (A / m^2)",
        legend = false
    )

    savefig(p, "graphs/current_0_special.pdf")
end
exit(1)
###########
# Current Graphs
i = 0
for curr_data in eachrow(current_data)
    global i
    println("[I] Current $(i) / $(size(curr_data, 1))")
    time_slice = curr_data[1]   # ms
    curr_data = curr_data[2:end]

    # Unnormalised
    p = plot(
        radius_vals, curr_data, 
        title = "Time $(time_slice) (ms)", 
        xlabel = "radius (m)", ylabel = "Current Density (A / m^2)",
        legend = false
    )

    savefig(p, "graphs/current_unnormalised_$(i).pdf")

    # Normalised
    M = maximum(abs.(curr_data))
    p = plot(
        radius_vals, curr_data / M, 
        title = "Time $(time_slice) (ms)", 
        xlabel = "radius (m)", ylabel = "Current Density (normalised)",
        legend = false
    )

    savefig(p, "graphs/current_normalised_$(i).pdf")

    i += 1
end

# Pressure Data
pressure_df = DataFrame(CSV.File("../isttok_p_profile.csv"))
pressure_data = Matrix(pressure_df)

# Pressure Graphs
i = 0
for pres_data in eachrow(pressure_data)
    global i
    println("[P] Pressure $(i) / $(size(pressure_data, 1))")
    time_slice = pres_data[1]
    pres_data = pres_data[2:end]

    # Unnormalised
    p = plot(
        radius_vals, pres_data, 
        title = "Time $(time_slice) (ms)", 
        xlabel = "radius (mm)", ylabel = "Pressure Density",
        legend = false
    )

    savefig(p, "graphs/pressure_unnormalised_$(i).pdf")

    # Normalised
    M = maximum(abs.(pres_data))
    p = plot(
        radius_vals, pres_data / M, 
        title = "Time $(time_slice) (ms)", 
        xlabel = "radius (mm)", ylabel = "Pressure Density (normalised)",
        legend = false
    )

    savefig(p, "graphs/pressure_normalised_$(i).pdf")

    i += 1
end