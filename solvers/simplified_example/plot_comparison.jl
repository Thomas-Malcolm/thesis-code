using Plots

a1,a2,a3 = [10.304680284464693, -2446.594091833618, 12.129380574926598]

f(a1,a2,a3,x) = a1 * x^2 + (a2 * x) / (a2 - 42.0) + a3


x_range = range(-200, 200, length = 1000)

p = plot(x_range, [ f(a1,a2,a3,xx) for xx in x_range ] , linecolor = :blue)

plot!(x_range, [ f(10.3, 1337.0, 12.0, xx) for xx in x_range], linecolor = :red )

savefig(p, "comparison.pdf")