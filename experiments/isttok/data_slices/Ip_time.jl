using CSV
using DataFrames

using Plots

cdata = Matrix(DataFrame(CSV.File("../isttok_ip_data.csv")))

time_slices = cdata[:,1]
jp = cdata[:,2]

p = plot(time_slices, jp, title = "Ip ISTTOK Data", legend = false)
xlabel!("Time (ms)")
ylabel!("Ip Total")

savefig(p, "graphs/jp_time_data.png")

