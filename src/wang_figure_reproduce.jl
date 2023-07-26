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

println("Current")
print("[")
for x in x_vals
    print("($(x), $(current_profile(x))), ")
end
print("]")

println()
println("current 2")

for x in x_vals
    print("($(x), $(current_vals(x))), ")
end
exit(1)


# x1 = x_range[20]
# x2 = x_range[100]
# x3 = last(x_range)
# x4 = x_range[150]
# x5 = x_range[70]

# c1 = current_vals(x1)
# c2 = current_vals(x2)
# c3 = current_vals(x3)
# c4 = current_vals(x4)
# c5 = current_vals(x5)

# println("Currents:")
# println("x $(x1) c $(c1)")
# println("x $(x2) c $(c2)")
# println("x $(x3) c $(c3)")
# println("x $(x4) c $(c4)")
# println("x $(x5) c $(c5)")


plot!(x_range, current_vals.(x_range))


β0 = minimum(2 * a1 * inner_iterate.(x_range, 0))

pressure_profile(x) = β0 - 2 * a1 * inner_iterate(x, 0)

maxPressure = maximum(abs.(pressure_profile.(x_range)))

pressure_vals(x) = pressure_profile(x) / maxPressure + 1
println()
println("Pressure")
print("[")
for x in x_vals
    print("($(x), $(pressure_vals(x))), ")
end
print("]")

plot!(x_range, pressure_vals.(x_range))



savefig("standard-a2.pdf")