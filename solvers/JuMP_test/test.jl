using JuMP
using NLopt
using Ipopt

using Plots

target_equation(x) = 10.3x^2 + 1337.0x / (1337.0 - 42.0) + 12.0
# target_equation(x) = 10.3x^2 + 1337.0x + 12.0


other_equation(a,b,c, x) = a * (x^2) + (( b ) / (b - 42.0)  ) * x + c
opto_eq(a,b,c, d) = abs(-d[2] + other_equation(a,b,c, d[1]))^2

# 4.0 ≤ x ≤ 6.0
d1 = (5.6, target_equation(5.6))
d2 = (4.9, target_equation(4.9))
d3 = (5.0, target_equation(5.0))
d4 = (4.3, target_equation(4.3))
d5 = (4.00001, target_equation(4.00001))
d6 = (5.4, target_equation(5.4))
d7 = (6.0, target_equation(6.0))

println(d1)
println(d2)
println(d3)

ds = [d1,d2,d3,d4,d5]


f(a,b,c) = sum(opto_eq(a,b,c,d) for d in ds)

function ∇f(G::AbstractVector{T}, a::T, b::T, c::T) where {T}


    G[1] = sum(
        2 * d[1]^2 * (
            -d[2] + other_equation(a,b,c, d[1])
        )
        for d in ds
    )

    G[2] = sum(
        2 * ( ( (-42.0) / (b - 42.0)^2 ) * d[1] ) * ( -d[2] + other_equation(a,b,c, d[1]) )
        for d in ds
    )

    G[3] = sum(
        2 * ( -d[2] + other_equation(a,b,c, d[1]) )
        for d in ds
    )

end


model = Model(Ipopt.Optimizer)#Model(NLopt.Optimizer)
# set_optimizer_attribute(model, "algorithm", :LD_MMA)

register(model, :my_eq, 3, f, ∇f)

@variable(model, 0 ≤ a ≤ 20)
@variable(model, 0 ≤ b ≤ 2000)
@variable(model, 0 ≤ c ≤ 30)

@NLobjective(model, Min, my_eq(a,b,c))

print(model)

print("Optimising...")
res = JuMP.optimize!(model)
print(res)

println("Obj val:", objective_value(model))
println("a,b,c=", [value(a), value(b), value(c)])