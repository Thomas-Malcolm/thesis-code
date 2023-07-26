using Bessels
using Roots
using Plots

# a1 = -0.04531
# a2 = -1.0808
# α = 2.1683
a1, a2, α = 0.01, 3.1, 5.566
# a1,a2,α = [-0.04, -1.1, 2.1]
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

eun(x, n) = (1 / μns[n]) * ( 
    a1 * x^2 * (
        cns[n] * besselj(2, μns[n] * x) + bessely(2, μns[n] * x)
    )
    +
    a2 * (
        cns[n] * besselj0(μns[n] * x) + bessely0(μns[n] * x)
    )
)

auns = [ (eun(x0 + 1, n) - eun(x0 - 1, n)) for n in 1:N ]

edn(x, n) = -0.5 * x^2 * bessely0(μns[n] * x) * bessely(2, μns[n] * x) - 0.5 * cns[n]^2 * x^2 * besselj0(μns[n] * x) * besselj(2, μns[n] * x) + cns[n] * (
    0.5 * x^2 * (
        besselj0(μns[n] * x) - besselj(2, μns[n] * x)
    )
    * bessely0(μns[n] * x)
    -
    (1/μns[n]) * x * besselj0(μns[n] * x) * bessely1(μns[n] * x)
)

adns = [ edn(x0 + 1, n) - edn(x0 - 1, n) for n in 1:N ]


inner_iterate(x, z) = x * sum(
    sum(
        (((-1)^(l-1) * 2 * auns[n]) / (k * vls[l] * adns[n] * (α^2 - λnls[n,l])))
        *
        (cns[n] * besselj1(μns[n] * x) + bessely1(μns[n] * x))
        *
        cos(vls[l] * z)
        for l in 1:N
    ) for n in 1:N
)

mag_field = [
    inner_iterate(x,z) for z in z_range, x in x_range
]




contour(x_range, z_range, mag_field)

current_profile(x) = -a1 * x + a2/x + (1/x)*(α^2) * inner_iterate(x,0)

maxCurrent = maximum(abs.(current_profile.(x_range)))

current_vals(x) = current_profile(x) / maxCurrent

x_vals = Set([rand(x_range) for _ in 1:50])

x_vals = Vector([x_vals...])



# curr1Vals = current_vals.(x_vals)
# curr2Vals = -curr1Vals

# p = scatter(
#     x_vals,
#     curr1Vals,
#     linecolor = :red
# )

# scatter!(
#     x_vals,
#     curr2Vals,
#     linecolor = :green
# )

# savefig(p, "banana.pdf")


function smoothly_invert(curr::Vector{Float64}, h::Int = 50)
    a = Animation()
    M = maximum(abs.(curr))

    targetCurrent = -curr

    for i in 1:h
        println("$(i)/50")

        newCurrent::Vector{Float64} = Vector{Float64}()

        for j in 1:length(curr)
            if curr[j] > 0
                push!(
                    newCurrent,
                    curr[j] - (abs(targetCurrent[j] - curr[j]) / h) * i
                )
            else
                push!(
                    newCurrent,
                    curr[j] + (abs(targetCurrent[j] - curr[j]) / h) * i
                )
            end
        end

        p = scatter(
            x_vals, 
            newCurrent,
            title = "banana",
            xlabel = "x",
            ylabel = "Current Vals",
            xlims = (4.0, 6.0),
            ylims = (-M, M)
        )

        frame(a)
    end

    gif(a, "current-inversion.gif"; fps = 5)
end

smoothly_invert(current_vals.(x_vals))

