using Bessels
using Roots
using Plots
using JuMP
using Ipopt
using NLopt

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
                ( 2 * tmpVal(n) * (-1)^(l-1)) 
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

# Optimisation functions - least squares

# Figure 1 data
ds = [(5.939698492462312, 0.1229276392921308), (4.1708542713567835, -0.010090512749661724), (4.874371859296482, 0.28282298944169054), (5.547738693467337, 0.28215565299241535), (5.477386934673367, 0.2994897328542968), (4.633165829145729, 0.20204291602150706), (4.472361809045226, 0.13344559596113667), (5.567839195979899, 0.2765369764260688), (4.613065326633166, 0.1939921156362315), (5.678391959798995, 0.2405237858455616), (4.743718592964824, 0.24308894237246098), (4.080402010050252, -0.052125947549691715), (5.768844221105527, 0.20463237273314155), (4.64321608040201, 0.20600583962216737), (4.492462311557789, 0.14250636321672805), (4.994974874371859, 0.3090397321598747), (5.728643216080402, 0.2212975119582323), (4.964824120603015, 0.3034920083780438), (4.7537688442211055, 0.24651516156413433), (5.949748743718593, 0.11775090351087136), (4.351758793969849, 0.07687980744116377), (5.135678391959799, 0.32584853663693464), (4.582914572864322, 0.18161617468716326), (4.100502512562814, -0.0428806277570339), (5.1457286432160805, 0.32646630666345067), (5.959798994974874, 0.11255824950729555), (5.748743718592965, 0.21310619071628578), (5.366834170854271, 0.3188711684065497), (5.175879396984925, 0.327838280105055), (4.894472361809045, 0.2879288208002983), (5.346733668341709, 0.32129361327193817), (5.698492462311558, 0.23304844923752568), (5.608040201005025, 0.2644379724026341), (4.241206030150754, 0.023468126297322756), (4.180904522613066, -0.005341526318123249), (5.708542713567839, 0.2292031771082979), (5.909547738693467, 0.13831566642315932), (5.9798994974874375, 0.10214468242514795), (5.025125628140704, 0.31390592685501517), (4.824120603015075, 0.2688328416357674), (5.71859296482412, 0.22528616594625292), (5.618090452261306, 0.26123511225158436)]
# Figure 2 data
# ds = [(4.844221105527638, 2.13578746422183), (4.562814070351759, -0.662168003246996), (4.874371859296482, 2.375173176933203), (4.050251256281407, 0.3166481951593722), (5.155778894472362, 2.79681072602925), (5.276381909547739, 2.0139644254424915), (4.592964824120603, -0.3976467524605849), (4.613065326633166, -0.20694960364205295), (5.678391959798995, -0.3605055430116794), (5.055276381909548, 3.04830972072795), (5.457286432160804, 0.52543966086227), (4.663316582914573, 0.30575580077424314), (4.341708542713568, -1.426169862063976), (4.080402010050252, 0.07127803472214478), (5.467336683417085, 0.45337894273394364), (5.407035175879397, 0.9160882551991394), (4.542713567839196, -0.8215800946618115), (5.236180904522613, 2.318832715319693), (4.64321608040201, 0.09578600492539502), (5.768844221105527, -0.2931606764799555), (4.311557788944723, -1.354733853244173), (4.834170854271357, 2.0495932396272782), (5.507537688442211, 0.1907939690585076), (5.557788944723618, -0.07074402365701316), (4.552763819095477, -0.7437054548311431), (5.045226130653266, 3.0490443128324034), (5.658291457286432, -0.3449283733081264), (4.703517587939698, 0.7357952165089913), (4.482412060301508, -1.2022240330842937), (4.381909547738694, -1.4591268274763385), (5.015075376884422, 3.023534250274462), (5.185929648241206, 2.6422052381354635), (5.4874371859296485, 0.31667293832937427), (5.346733668341709, 1.425338773619112), (4.4120603015075375, -1.4337828040912477), (4.231155778894473, -0.9994574991406685), (5.648241206030151, -0.3324434822795202), (4.201005025125628, -0.8171328592320788), (5.8291457286432165, -0.15052516845479402), (5.5879396984924625, -0.18834297337570982), (5.226130653266332, 2.389284950904071), (5.9798994974874375, 0.3790723159534339), (5.71859296482412, -0.3563718152581148), (4.8040201005025125, 1.7741873353335769), (5.105527638190955, 2.9765232903559626)]


function residual(a1, a2, α, d)
    M = maxIn(a1, a2, α)

    d[2] - unIn(d[1], a1, a2, α)
end


sum_residuals(a1, a2, α) = sum(residual(a1, a2, α, d)^2 for d in ds)

function f(a1, a2, α)
    println("--new-iteration---")
    println("a1,a2,α= $(a1), $(a2), $(α)")

    # println("RESIDUALS")
    # # for d in ds
    # #     resid = residual(a1, a2, α, d)^2
    # #     println("  diff: $(resid)")
    # #     println("--")
    # # end
    println("total residual: ", sum(residual(a1,a2,α,d)^2 for d in ds))

    sum(residual(a1,a2,α,d)^2 for d in ds)
end


function ∇f(G::AbstractVector{T}, a1::T, a2::T, α::T) where T

    M = maxIn(a1, a2, α)

    G[1] = sum(
        -2 * (
            d[2] - unIn(d[1], a1, a2, α)
        ) * InDa1(d[1], a1, a2, α)
        for d in ds
    )

    G[2] = sum(
        -2 * (
            d[2] - unIn(d[1], a1, a2, α)
        ) * InDa2(d[1], a1, a2, α)
        for d in ds
    )

    G[3] = sum(
        -2* (
            d[2] - unIn(d[1], a1, a2, α)
        ) * InDalpha(d[1], a1, a2, α)
        for d in ds
    )

    println("Der wrt a1: $(G[1])")
    println("Der wrt a2: $(G[2])")
    println("Der wrt α : $(G[3])")
end


# Optimising
println("Optimising")

model = Model(NLopt.Optimizer)
# model = Model(Ipopt.Optimizer)


set_optimizer_attribute(model, "algorithm", :LD_MMA)
# local_optimizer = NLopt.Opt(:LD_LBFGS, 3)
# local_optimizer.xtol_abs = 0
# local_optimizer.ftol_abs = 1e-5
# set_optimizer_attribute(model, "local_optimizer", local_optimizer)
set_optimizer_attribute(model, "ftol_abs", 1e-5)
set_optimizer_attribute(model, "xtol_abs", 0)       # this seems to not be 0 despite docs saying it defaults to that?

register(model, :pls_work, 3, f, ∇f)

# Figure 1
@variable(model, -0.06 ≤ a1 ≤ -0.04)
@variable(model, -1.5 ≤ a2 ≤ -0.5)
@variable(model, 2.0 ≤ α ≤ 3.0)
# Figure 2
# @variable(model, 0.0 ≤ a1 ≤ 0.05)
# @variable(model, 2.8 ≤ a2 ≤ 3.2)
# @variable(model, 5.0 ≤ α ≤ 5.8)

# Figure 1
# [-0.04531, -1.0808, 2.1683]
set_start_value(a1, -0.0475)
set_start_value(a2, -1.06)
set_start_value(α, 2.167)

# Figure 2
# [0.01, 3.1, 5.566]
# set_start_value(a1, 0.02)
# set_start_value(a2, 3.0)
# set_start_value(α, 5.2)



@NLobjective(model, Min, pls_work(a1, a2, α))

print(model)

println("Optimising...")
JuMP.optimize!(model)
print(solution_summary(model))

println("Obj val:", objective_value(model))
print("Termination:", termination_status(model))
println("a1,a2,α=", [value(a1), value(a2), value(α)])

a1,a2,α = [value(a1), value(a2), value(α)]


println("residuals:")
for d in ds
    println("- $(residual(a1, a2, α, d)^2)")
end




## SEEMS THAT THIS WORKS JUST AS WELL - MAYBE MOVE TO THIS METHOD, SEEING AS THE CURRENT
##  WE WILL RECEIVE WILL OF COURSE NOT BE NORMALISED

# Graph Making


## Magnetic Field
mag_field = [
    ψ(a1, a2, α, x, z) for z in z_range, x in x_range
]
contour(x_range, z_range, mag_field)

## Current
curr_vals = [In(xx, a1, a2, α) for xx in x_range]
plot!(x_range, curr_vals)

## Pressure
β0 = minimum(2 * a1 * ψ.(a1, a2, α, x_range, 0))

pressure_profile(x) = β0 - 2 * a1 * ψ.(a1, a2, α, x, 0)

maxPressure = maximum(abs.(pressure_profile.(x_range)))

pressure_vals(x) = pressure_profile(x) / maxPressure + 1

plot!(x_range, pressure_vals.(x_range))
savefig("figure-1-no-local-test.pdf")
