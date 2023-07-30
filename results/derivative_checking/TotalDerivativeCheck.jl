include("../../src/ACSimulation.jl")
using .ACSimulation
using Plots

# Default
cfg = Config()
x_range = cfg.x_range
z_range = cfg.z_range

# Figure 1
params = Parameters(-0.04531, -1.0808, 2.1683)
ds = Data([(4.562814070351759, 0.5271319912938993), (4.1708542713567835, -0.03071488391231543), (4.442211055276382, 0.3641464943325205), (4.874371859296482, 0.8608953285081333), (4.050251256281407, -0.200707612808032), (5.849246231155779, 0.5116083973086332), (4.693467336683417, 0.6853595225343578), (6.0, 0.2792101836871248), (4.633165829145729, 0.6150058837311638), (5.276381909547739, 0.9954062001816789), (4.733668341708543, 0.729349784351659), (4.743718592964824, 0.7399473971813753), (5.899497487437186, 0.43644828097097393), (4.080402010050252, -0.15866809423157385), (5.35678391959799, 0.9744431765531748), (5.467336683417085, 0.9182343507501194), (4.391959798994975, 0.29246728959561613), (4.0703517587939695, -0.17269812613700056), (4.311557788944723, 0.17497300919682354), (4.14070351758794, -0.07377749046144401), (4.954773869346734, 0.9177227304059346), (4.834170854271357, 0.8272421593617327), (5.557788944723618, 0.8504227415781456), (5.035175879396985, 0.9599844963148867), (5.165829145728643, 0.9967729184103393), (5.738693467336684, 0.6612569751687866), (5.809045226130653, 0.5687994533183357), (4.703517587939698, 0.6965887940507236), (5.819095477386934, 0.5547786010884644), (5.21608040201005, 1.0), (4.582914572864322, 0.5528281724139369), (5.517587939698492, 0.8828491946944586), (4.321608040201005, 0.18976790719521844), (5.376884422110553, 0.9665403309747135), (5.256281407035176, 0.9979844333857804), (5.1457286432160805, 0.9937428314320712), (5.778894472361809, 0.6096737484118069), (5.1959798994974875, 0.9994640685431111), (4.4623115577889445, 0.3922702205715146), (4.090452261306533, -0.14461259382704408), (4.231155778894473, 0.056720276514437375), (5.206030150753769, 0.9998591389659488), (5.869346733668341, 0.48197986791545677), (4.221105527638191, 0.04203981316528289)])

# Figure 2
# params = Parameters(0.01, 3.1, 5.566)
# need data

## Testing InDa1
println("A1 partial")
Δh = 0.00001

a1_range = range(-0.06, -0.02, length = 100)
calcDer = [ residualDa1(cfg, Parameters(a1, params.a2, params.α), ds) for a1 in a1_range ]

manualDer = [ ( residuals(cfg, Parameters(a1 + Δh, params.a2, params.α), ds) - residuals(cfg, Parameters(a1, params.a2, params.α), ds) ) / Δh for a1 in a1_range ]
residualVals = [ residuals(cfg, Parameters(a1, params.a2, params.α), ds) for a1 in a1_range ]

p = plot(a1_range, calcDer, title = "Partial wrt a1", label = "Analytic")
plot!(a1_range, manualDer, label = "Manual")
plot!(a1_range, residualVals, label = "Residual")
vline!([params.a1], label = "Expected")
savefig(p, "partial-a1-residuals-fig-1.pdf")



# ## Testing InDa2
println("A2 partial")
Δh = 0.00001

a2_range = range(-1.2, -0.8, length = 100)
calcDer = [ residualDa2(cfg, Parameters(params.a1, a2, params.α), ds) for a2 in a2_range ]

manualDer = [ ( residuals(cfg, Parameters(params.a1, a2 + Δh, params.α), ds) - residuals(cfg, Parameters(params.a1, a2, params.α), ds) ) / Δh for a2 in a2_range ]
residualVals = [ residuals(cfg, Parameters(params.a1, a2, params.α), ds) for a2 in a2_range ]

p = plot(a2_range, calcDer, title = "Partial wrt a2", label = "Analytic")
plot!(a2_range, manualDer, label = "Manual")
plot!(a2_range, residualVals, label = "Residual")
vline!([params.a2], label = "Expected")
savefig(p, "partial-a2-residuals-fig-1.pdf")


## Testing InDalpha
println("α partial")
Δh = 0.00001

a3_range = range(1.9, 2.2, length = 100)
calcDer = [ residualDalpha(cfg, Parameters(params.a1, params.a2, α), ds) for α in a3_range ]

manualDer = [ ( residuals(cfg, Parameters(params.a1, params.a2, α + Δh), ds) - residuals(cfg, Parameters(params.a1, params.a2, α), ds) ) / Δh for α in a3_range ]
residualVals = [ residuals(cfg, Parameters(params.a1, params.a2, α), ds) for α in a3_range ]

p = plot(a3_range, calcDer, title = "Partial wrt alpha", label = "Analytic")
plot!(a3_range, manualDer, label = "Manual")
plot!(a3_range, residualVals, label = "Residual")
vline!([params.α], label = "Expected")
savefig(p, "partial-alpha-residuals-fig-1.pdf")