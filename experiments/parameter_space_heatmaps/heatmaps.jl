using Bessels
using Roots
using PlotlyJS
using DataFrames

x0 = 5.0
ρ = 1.0
k = ρ

h = 200
N = 8

x_range = range(x0 - ρ, x0 + ρ, length = h)
z_range = range(-ρ, ρ, length = h)


vls = [
    (l + 0.5) * π / k for l in 0:N-1
]

zero_function(μn) = -besselj1(μn * (x0 + 1)) * ( bessely1(μn * (x0 - 1)) / besselj1(μn * (x0 - 1)) ) + bessely1(μn * (x0 + 1))

zeros = find_zeros(
    zero_function, (0, 2π*N)
)

true_zeros = findall(
    μn -> abs(zero_function(μn)) < 10e-10, zeros
)

μns = zeros[splice!(true_zeros, 1:N)]

cn(μn) = -bessely1(μn * (x0 - 1)) / besselj1(μn * (x0 - 1))

cns = cn.(μns)

λnls = [ vls[l]^2 + μns[n]^2 for n in 1:N, l in 1:N ]







# Things after this are dependent on a1,a2 and/or α

eun(a1, a2, x, n) = (1 / μns[n]) * ( 
    a1 * x^2 * (
        cns[n] * besselj(2, μns[n] * x) + bessely(2, μns[n] * x)
    )
    +
    a2 * (
        cns[n] * besselj0(μns[n] * x) + bessely0(μns[n] * x)
    )
)

edn(x, n) = -0.5 * x^2 * bessely0(μns[n] * x) * bessely(2, μns[n] * x) - 0.5 * cns[n]^2 * x^2 * besselj0(μns[n] * x) * besselj(2, μns[n] * x) + cns[n] * (
    0.5 * x^2 * (
        besselj0(μns[n] * x) - besselj(2, μns[n] * x)
    )
    * bessely0(μns[n] * x)
    -
    (1/μns[n]) * x * besselj0(μns[n] * x) * bessely1(μns[n] * x)
)


function ψ(
    a1,
    a2,
    α,
    x,
    z
)

    auns = [ (eun(a1, a2, x0 + 1, n) - eun(a1, a2, x0 - 1, n)) for n in 1:N ]
    adns = [ (edn(x0 + 1, n) - edn(x0 - 1, n)) for n in 1:N ]

    x * sum(
        sum(
            (((-1)^(l-1) * 2 * auns[n]) / (k * vls[l] * adns[n] * (α^2 - λnls[n,l])))
            *
            (cns[n] * besselj1(μns[n] * x) + bessely1(μns[n] * x))
            *
            cos(vls[l] * z)
            for l in 1:N
        ) for n in 1:N
    )
end

function maxIn(a1, a2, α)
    current_profile(x) = -a1 * x + a2/x + (1/x) * (α^2) * ψ(a1, a2, α, x, 0)
    maximum(abs.(current_profile.(x_range)))
end

function unIn(xx, a1, a2, α)
    -a1 * xx + a2/xx + (1/xx) * (α^2) * ψ(a1, a2, α, xx, 0)
end


# Figure 1 data
current_data = [(4.874371859296482, 0.8608953285081333), (4.050251256281407, -0.200707612808032), (5.2462311557788945, 0.9988786605043959), (4.633165829145729, 0.6150058837311638), (4.733668341708543, 0.729349784351659), (4.743718592964824, 0.7399473971813753), (4.542713567839196, 0.5009892797476437), (5.628140703517588, 0.785217595998286), (5.417085427135678, 0.9476015319865697), ]#(5.0954773869346734, 0.9819407200465001), (4.924623115577889, 0.8980929903036065), (5.336683417085427, 0.9812856398146681), (4.311557788944723, 0.17497300919682354), (5.597989949748744, 0.8144653808836534), (5.306532663316583, 0.989547707387002), (4.994974874371859, 0.9406974385817826), (4.432160804020101, 0.34995927003114), (5.035175879396985, 0.9599844963148867), (4.8542713567839195, 0.8444898451426189), (5.809045226130653, 0.5687994533183357), (5.738693467336684, 0.6612569751687866), (4.703517587939698, 0.6965887940507236), (5.949748743718593, 0.35842631803105507), (4.351758793969849, 0.23401728131562705), (5.819095477386934, 0.5547786010884644), (4.5125628140703515, 0.46097149860085346), (5.075376884422111, 0.975557797775825), (4.522613065326633, 0.4744153415301792), (4.271356783919598, 0.11574061789509066), (5.919597989949748, 0.40549665777593163), (4.422110552763819, 0.335693459085735), (5.125628140703518, 0.9897405057750671), (4.241206030150754, 0.07143549517699091), (5.708542713567839, 0.6976800041654228), (4.984924623115578, 0.935298410971659), (5.9798994974874375, 0.31092196608679995), (5.63819095477387, 0.7750356576165367), (4.944723618090452, 0.9114066188298575), (5.5376884422110555, 0.8670833594503413), (5.71859296482412, 0.6857568694238816), (4.8040201005025125, 0.7998473179072897), (5.105527638190955, 0.9847783694102789)]
# pressure_data = [(4.874371859296482, 0.8842180169916798), (4.050251256281407, 0.045047022603066145), (5.2462311557788945, 0.9944947697510185)]#, (4.633165829145729, 0.6677242353053197), (4.733668341708543, 0.7671345281959976), (4.743718592964824, 0.7764750744620412), (4.542713567839196, 0.5712834796545643), (5.628140703517588, 0.701703169499505), (5.417085427135678, 0.9192112296402145),]# (5.0954773869346734, 0.9899015111953064), (4.924623115577889, 0.9175198007193006), (5.336683417085427, 0.9660180899452057), (4.311557788944723, 0.3109807737447373), (5.597989949748744, 0.7409298258787953), (5.306532663316583, 0.9783100488351647), (4.994974874371859, 0.9552768912930163), (4.432160804020101, 0.4477805457156936), (5.035175879396985, 0.9719387358848759), (4.8542713567839195, 0.8695197222669107), (5.809045226130653, 0.40750941238711424), (5.738693467336684, 0.5340973985226456), (4.703517587939698, 0.7383896841785181), (5.949748743718593, 0.11345663741481637), (4.351758793969849, 0.35634328528816683), (5.819095477386934, 0.38818557228569806), (4.5125628140703515, 0.538086593879546), (5.075376884422111, 0.9848732519759199), (4.522613065326633, 0.5492011090859339), (4.271356783919598, 0.2663530573224481), (5.919597989949748, 0.1800864434396765), (4.422110552763819, 0.4363680021231472), (5.125628140703518, 0.9955903451610858), (4.241206030150754, 0.23358350581812226), (5.708542713567839, 0.5835857931461892), (4.984924623115578, 0.950543213944616), (5.9798994974874375, 0.045632719089617524), (5.63819095477387, 0.6880211276370274), (4.944723618090452, 0.9293934020516461), (5.5376884422110555, 0.8112764379664622), (5.71859296482412, 0.5674085050307947), (4.8040201005025125, 0.8296156504885825), (5.105527638190955, 0.9920476165029324)]

function residual_In(a1, a2, α, d)
    M = maxIn(a1, a2, α)

    d[2] - (1/M) * unIn(d[1], a1, a2, α)
end

# function residual_pressure(a1, a2, α, d)
#     β0 = minimum([2 * a1 * ψ(a1, a2, α, x, 0) for x in x_range])

#     pressure_profile(x) = β0 - 2 * a1 * ψ(a1, a2, α, x, 0)

#     maxPressure = maximum(abs.(pressure_profile.(x_range)))

#     pressure_vals(x) = pressure_profile(x) / maxPressure + 1

#     d[2] - pressure_vals(d[1])
# end

sum_residuals(a1, a2, α) = sum(residual_In(a1, a2, α, d)^2 for d in current_data)# + 
    # sum(residual_pressure(a1, a2, α, d)^2 for d in pressure_data)


# Figure 1 approximation starting values

# Coarse
# x1_vals = range(-0.06, 0.06, length = 40)
# x2_vals = range(-2, 2, length = 40)
# x3_vals = range(1.5, 2.5, length = 40)

# Specific
# x1_vals = range(-0.046, -0.044, length = 40)
# x2_vals = range(-1.1, -1.0, length = 40)
# x3_vals = range(2.16, 2.17, length = 40)

x1_vals = range(-0.048, -0.045, length = 40)
x2_vals = range(-1.1, -1.0, length = 40)
x3_vals = range(2.16, 2.17, length = 40)


# Varying a1 and a2, with a fixed α
println("a1, a2, fixed α")

heat_vals = [
    sum_residuals(a1v, a2v, 2.1683)
    for a2v in x2_vals, a1v in x1_vals
]

df1 = DataFrame(
    x = [-0.04531 for _ in 1:40],
    y = x2_vals
)

df2 = DataFrame(
    x = x1_vals,
    y =  [-1.0808 for _ in 1:40]
)

hm1 = plot([
        heatmap(
            x = x1_vals,
            y = x2_vals,
            z = heat_vals,
        ),
        scatter(
            x = df1.x,
            y = df1.y,
        ),
        scatter(
            x = df2.x,
            y = df2.y,
        )
    ],
    Layout(
        xaxis_title = "a1",
        yaxis_title = "a2",
        showlegend = false
    )
)

savefig(hm1, "heat-1-blue-only.pdf")
exit(1)

# Varying a1 and α, with a fixed a2
println("a1, α, fixed a2")
heat_vals = [
    sum_residuals(a1v, -1.0808, αv)
    for αv in x3_vals, a1v in x1_vals
]

df1 = DataFrame(
    x = [-0.04531 for _ in 1:40],
    y = x3_vals
)

df2 = DataFrame(
    x = x1_vals,
    y = [2.1683 for _ in 1:40]
)


hm2 = plot([
        heatmap(
            x = x1_vals,
            y = x3_vals,
            z = heat_vals,
        ),
        scatter(x = df1.x, y = df1.y),
        scatter(x = df2.x, y = df2.y)
    ],
    Layout(
        xaxis_title = "a1",
        yaxis_title = "alpha",
        showlegend = false
    )
)

savefig(hm2, "heat-2-entire.pdf")

# Varying a2 and α, with a fixed a1
println("a2, α, fixed a1")
heat_vals = [
    sum_residuals(-0.04531, a2v, αv)
    for αv in x3_vals, a2v in x2_vals
]

df1 = DataFrame(
    x = [-1.0808 for _ in 1:40],
    y = x3_vals
)

df2 = DataFrame(
    x = x2_vals,
    y =  [2.1683 for _ in 1:40]
)

# hm3 = plot(heatmap(
#     x = x2_vals,
#     y = x3_vals,
#     z = heat_vals
# ))

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

savefig(hm3, "heat-3-entire.pdf")
# p = [hm1 hm2 ; hm3 nothing]
# relayout!(p)
# p = plot(hm1, hm2, hm3, layout = @layout [a b ; c d])
# p = make_subplots(
#     rows = 2, cols = 2,
#     specs = [
#         Spec() Spec() ; Spec()  missing
#     ]
# )
# add_trace!(p, hm1, row = 1, col = 1)
# add_trace!(p, hm2, row = 1, col = 2)
# add_trace!(p, hm3, row = 2, col = 1)

# savefig(p, "conjoined.pdf")