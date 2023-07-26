using NLsolve
using Optim
using LineSearches

using Plots

target_equation(x) = 10.3x^2 + 1337.0x / (1337.0 - 42.0) + 12.0
# target_equation(x) = 10.3x^2 + 1337.0x + 12.0



#opto_eq(Z, d) = abs(-d[2] + (Z[1] * (d[1]^2) + ((Z[2] / (Z[2] - 42.0)) * d[1]) + Z[3]))^2

# opto_eq(Z, d) = abs(-d[2] +  Z[1] * (d[1]^2) + (  (Z[2]) / (Z[2] - 42.0)  ) * d[1] + Z[3])^2
other_equation(Z, x) = Z[1] * (x^2) + (  (Z[2]) / (Z[2] - 42.0)  ) * x + Z[3]
opto_eq(Z, d) = abs(-d[2] + other_equation(Z, d[1]))^2

# 4.0 ≤ x ≤ 6.0
d1 = (5.6, target_equation(5.6))
d2 = (4.9, target_equation(4.9))
d3 = (5.0, target_equation(5.0))

println(d1)
println(d2)
println(d3)

# function f!(F, x)

#     F[1] = opto_eq(x, d1)
#     F[2] = opto_eq(x, d2)
#     F[3] = opto_eq(x, d3)
# end

# a = Animation()
# x_range = range(4.0, 6.0, length = 200)

# expected_vals = target_equation.(x_range)

# counter = 0

function to_minimise(Z)

    # global counter

    objVal = opto_eq(Z, d1) + opto_eq(Z, d2) + opto_eq(Z, d3)

    # yVals = [other_equation(Z, xx) for xx in x_range]

    # p = plot(x_range, yVals, title="($(counter)) a: $(round(Z[1]; digits = 3)) b: $(round(Z[2]; digits = 3)) c: $(round(Z[3]; digits = 3))", ylim = (0, 700))
    # plot!(x_range, expected_vals, linecolor = :red)
    # frame(a)

    # counter += 1

    return objVal
end



# TODO: try different linesearch algorithms. The asymptotic nature of it seems to be messing with the step sizes.


function g!(G, Z)

    G[1] = 2 * d1[1]^2 * ( -d1[2] + other_equation(Z, d1[1]) ) +
            2 * d2[1]^2 * ( -d2[2] + other_equation(Z, d2[1]) ) +
            2 * d3[1]^2 * ( -d3[2] + other_equation(Z, d3[1]) )

    G[2] = 2 * ( ( (-42.0) / (Z[2] - 42.0)^2 ) * d1[1] ) * ( -d1[2] + other_equation(Z, d1[1]) ) +
            2 * ( ( (-42.0) / (Z[2] - 42.0)^2 ) * d2[1] ) * ( -d2[2] + other_equation(Z, d2[1]) ) +
            2 * ( ( (-42.0) / (Z[2] - 42.0)^2 ) * d3[1] ) * ( -d3[2] + other_equation(Z, d3[1]) )

    G[3] = 2 * ( -d1[2] + other_equation(Z, d1[1]) ) +
            2 * ( -d2[2] + other_equation(Z, d2[1]) ) +
            2 * ( -d3[2] + other_equation(Z, d3[1]) )
end

# to_minimise(Z) = opto_eq(Z, d1) + opto_eq(Z, d2) + opto_eq(Z, d3)

# Minimise it
# res = optimize(to_minimise, [5.0, 2.0, 30.0], BFGS(), Optim.Options(
#     iterations = 5000,
#     #store_trace = true,
#     #show_trace = true
# ); autodiff = :forward)


# Should we specify nonlinear preconditions?
# res = optimize(
#     to_minimise,
#     g!,
#     [5.0, 2.0, 30.0], 
#     NGMRES(
#         ; linesearch = LineSearches.HagerZhang()
#     ),
#     Optim.Options(
#         iterations = 10000
#     )
# )
println("Optimising")
res = optimize(
    to_minimise,
    g!,
    [5.0, 2.0, 30.0], 
    Newton(
        ; 
        alphaguess = LineSearches.InitialStatic(),
        linesearch = LineSearches.BackTracking()#linesearch = LineSearches.MoreThuente()
    ),
    Optim.Options(
        iterations = 10000
    )
)

# gif(a, "newton-back_tracking.gif")

println(res)
println("a,b,c:")
println(Optim.minimizer(res))