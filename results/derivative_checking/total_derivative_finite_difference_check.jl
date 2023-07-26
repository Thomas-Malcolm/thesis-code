using Bessels
using Roots
using Plots
# using PlotlyJS
using Optim
using LineSearches
# using JuMP
# using Ipopt
# using NLopt

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
    -a1 * xx + a2/xx + (1/xx)*(α^2) * ψ(a1, a2, α, xx, 0)
end

function In(xx, a1, a2, α)
    current_profile(x) = -a1 * x + a2/x + (1/x) * (α^2) * ψ(a1, a2, α, x, 0)
    maxCurrent = maximum(abs.(current_profile.(x_range)))

    current_profile(xx) / maxCurrent
end

# Partial derivative of current wrt a1
function InDa1(xx, a1, a2, α)

    tmpVal(n) = (1/μns[n]) * (
        (x0+1)^2 * (
            cns[n] * besselj(2, μns[n] * (x0 + 1)) + bessely(2, μns[n] * (x0 + 1))
        ) 
        -
        (x0-1)^2 * (
            cns[n] * besselj(2, μns[n] * (x0 - 1)) + bessely(2, μns[n] * (x0 - 1))
        )
    )

    adns = [ (edn(x0 + 1, n) - edn(x0 - 1, n)) for n in 1:N ]

    -xx + α^2 * sum(
        sum(
            (
                ( (-1)^(l-1) * 2 * tmpVal(n)) 
                / 
                ( k * vls[l] * adns[n] * (α^2 - λnls[n,l]) )
            )
            *
            (cns[n] * besselj1(μns[n] * xx) + bessely1(μns[n] * xx) )
            for l in 1:N
        ) for n in 1:N
    )
end


# Partial derivative of current wrt a2
function InDa2(xx, a1, a2, α)

    tmpVal(n) = (1/μns[n]) * (
        (
            cns[n] * besselj0(μns[n] * (x0 + 1)) + bessely0(μns[n] * (x0 + 1))
        ) 
        -
        (
            cns[n] * besselj0(μns[n] * (x0 - 1)) + bessely0(μns[n] * (x0 - 1))
        )
    )

    adns = [ (edn(x0 + 1, n) - edn(x0 - 1, n)) for n in 1:N ]

    (1/xx) + α^2 * sum(
        sum(
            (
                ( (-1)^(l-1) * 2 * tmpVal(n)) 
                / 
                ( k * vls[l] * adns[n] * (α^2 - λnls[n,l]) )
            )
            *
            (cns[n] * besselj1(μns[n] * xx) + bessely1(μns[n] * xx) )
            for l in 1:N
        ) for n in 1:N
    )
end

# Partial derivative of current wrt α
function InDalpha(xx, a1, a2, α)

    auns = [ (eun(a1, a2, x0 + 1, n) - eun(a1, a2, x0 - 1, n)) for n in 1:N ]
    adns = [ (edn(x0 + 1, n) - edn(x0 - 1, n)) for n in 1:N ]

    part1 = 2 * α * (1/xx) * ψ(a1,a2,α,xx,0)

    part2 = α^2 * sum(
        sum(
            (
                ((-1)^(l) * 4 * auns[n] * α)
                / 
                (k * vls[l] * adns[n] * (α^2 - λnls[n,l])^2)
            )
            *
            (cns[n] * besselj1(μns[n] * xx) + bessely1(μns[n] * xx))
            for l in 1:N
        ) for n in 1:N
    )

    return part1 + part2
end

# Figure 1 data
ds = [(4.562814070351759, 0.5271319912938993), (4.1708542713567835, -0.03071488391231543), (4.442211055276382, 0.3641464943325205), (4.874371859296482, 0.8608953285081333), (4.050251256281407, -0.200707612808032), (5.849246231155779, 0.5116083973086332), (4.693467336683417, 0.6853595225343578), (6.0, 0.2792101836871248), (4.633165829145729, 0.6150058837311638), (5.276381909547739, 0.9954062001816789), (4.733668341708543, 0.729349784351659), (4.743718592964824, 0.7399473971813753), (5.899497487437186, 0.43644828097097393), (4.080402010050252, -0.15866809423157385), (5.35678391959799, 0.9744431765531748), (5.467336683417085, 0.9182343507501194), (4.391959798994975, 0.29246728959561613), (4.0703517587939695, -0.17269812613700056), (4.311557788944723, 0.17497300919682354), (4.14070351758794, -0.07377749046144401), (4.954773869346734, 0.9177227304059346), (4.834170854271357, 0.8272421593617327), (5.557788944723618, 0.8504227415781456), (5.035175879396985, 0.9599844963148867), (5.165829145728643, 0.9967729184103393), (5.738693467336684, 0.6612569751687866), (5.809045226130653, 0.5687994533183357), (4.703517587939698, 0.6965887940507236), (5.819095477386934, 0.5547786010884644), (5.21608040201005, 1.0), (4.582914572864322, 0.5528281724139369), (5.517587939698492, 0.8828491946944586), (4.321608040201005, 0.18976790719521844), (5.376884422110553, 0.9665403309747135), (5.256281407035176, 0.9979844333857804), (5.1457286432160805, 0.9937428314320712), (5.778894472361809, 0.6096737484118069), (5.1959798994974875, 0.9994640685431111), (4.4623115577889445, 0.3922702205715146), (4.090452261306533, -0.14461259382704408), (4.231155778894473, 0.056720276514437375), (5.206030150753769, 0.9998591389659488), (5.869346733668341, 0.48197986791545677), (4.221105527638191, 0.04203981316528289)]


function residual(a1, a2, α, d)
    M = maxIn(a1, a2, α)

    d[2] - (1/M) * unIn(d[1], a1, a2, α)
end


function f(a1, a2, α)
    sum(residual(a1,a2,α,d)^2 for d in ds)
end

function totalDerA1(a1, a2, α)
    M = maxIn(a1, a2, α)

    sum(
        (-2 / M) * (
            d[2] - (1/M) * unIn(d[1], a1, a2, α)
        ) * InDa1(d[1], a1, a2, α)
        for d in ds
    )
end

## Testing InDa1
println("A1 partial")
Δh = 0.000001
a1_range = range(-0.06, -0.02, length = 50)
calcDer = [ totalDerA1(a1, -1.0808, 2.1683 ) for a1 in a1_range ]
manualDer = [ (f(a1 + Δh, -1.0808, 2.1683) - f(a1, -1.0808, 2.1683) ) / Δh for a1 in a1_range ]
residualVals = [ f(a1, -1.0808, 2.1683) for a1 in a1_range ]

p = plot(a1_range, calcDer, title = "Partial wrt a1", label = "Analytic")
plot!(a1_range, manualDer, label = "Manual")
plot!(a1_range, residualVals, label = "Residual")
vline!([-0.04531], label = "Expected")
savefig(p, "partial-a1-residuals.pdf")

function totalDerA2(a1, a2, α)

    M = maxIn(a1, a2, α)

    sum(
        (-2 / M) * (
            d[2] - (1/M) * unIn(d[1], a1, a2, α)
        ) * InDa2(d[1], a1, a2, α)
        for d in ds
    )
end

## Testing InDa2
println("A2 partial")
a2_range = range(-1.2, -0.8, length = 50)
calcDer = [ totalDerA2(-0.04531, a2, 2.1683 ) for a2 in a2_range ]
manualDer = [ (f(-0.04531, a2 + Δh, 2.1683) - f(-0.04531, a2, 2.1683) ) / Δh for a2 in a2_range ]
residualVals = [ f(-0.04531, a2, 2.1683) for a2 in a2_range ]

p = plot(a2_range, calcDer, title = "Partial wrt a2", label = "Analytic")
plot!(a2_range, manualDer, label = "Manual")
plot!(a2_range, residualVals, label = "Residual")
vline!([-1.0808], label = "Expected")
savefig(p, "partial-a2-residuals.pdf")


function totalDerAlpha(a1, a2, α)

    M = maxIn(a1, a2, α)

    sum(
        (-2 / M) * (
            d[2] - (1/M) * unIn(d[1], a1, a2, α)
        ) * InDalpha(d[1], a1, a2, α)
        for d in ds
    )
end

println("α partial")
a3_range = range(1.9, 2.2, length = 50)
calcDer = [ totalDerAlpha(-0.04531, -1.0808, alpha ) for alpha in a3_range ]
manualDer = [ (f(-0.04531, -1.0808, alpha + Δh) - f(-0.04531, -1.0808, alpha) ) / Δh for alpha in a3_range ]
residualVals = [ f(-0.04531, -1.0808, alpha) for alpha in a3_range ]

p = plot(a3_range, calcDer, title = "Partial wrt α", label = "Analytic")
plot!(a3_range, manualDer, label = "Manual")
plot!(a3_range, residualVals, label = "Residual")
vline!([2.1683], label = "Expected")
savefig(p, "partial-alpha-residuals.pdf")

