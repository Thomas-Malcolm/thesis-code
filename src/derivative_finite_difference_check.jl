using Bessels
using Roots
# using Plots
using PlotlyJS
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


## Testing InDa1
somea1 = (-0.045, 0.05, 0.1, -0.06)
calcDer = [ InDa1(5.0, a1, -1.0808, 2.1683 ) for a1 in somea1 ]
manualDer = [ (unIn(5.0, a1 + 0.0001, -1.0808, 2.1683) - unIn(5.0, a1, -1.0808, 2.1683) ) / (0.0001) for a1 in somea1 ]

println("InDa1")
println("calculated derivatives:")
println(calcDer)
println("finite difference derivatives:")
println(manualDer)




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


## Testing InDa2
somea2 = (-1.08, -1.2, -1.0, 2.0, 1.73)
calcDer = [ InDa2(5.0, -0.04531, a2, 2.1683 ) for a2 in somea2 ]
manualDer = [ (unIn(5.0, -0.04531, a2 + 0.0001, 2.1683) - unIn(5.0, -0.04531, a2, 2.1683) ) / (0.0001) for a2 in somea2 ]

println("InDa2")
println("calculated derivatives:")
println(calcDer)
println("finite difference derivatives:")
println(manualDer)



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


## Testing InDa2
somea3 =  range(2.1683 - 0.005, 2.1683 + 0.005, length = 10) #[2.1683 + 0.0001 * 1:10]
calcDer = [ InDalpha(5.0, -0.04531, -1.0808, alpha ) for alpha in somea3 ]
manualDer = [ (unIn(5.0, -0.04531, -1.0808, alpha + 0.0001) - unIn(5.0, -0.04531, -1.0808, alpha) ) / (0.0001) for alpha in somea3 ]

println("InDalpha")
println("calculated derivatives:")
println(calcDer)
println("finite difference derivatives:")
println(manualDer)
