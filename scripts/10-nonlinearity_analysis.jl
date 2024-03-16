#=
Цель программы: проанализировать нелинейность f(x,a,b,c) при конкретных значениях a,b,c.
Под этим подразумеваются точки и значения экстремумов.

Значения a,b,c: a = 0.32, b = 0.79, c = 1.166
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie
using Zygote
using Optim

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001

from_c_offset = 0.06
x_start, x_end, x_N = -c-from_c_offset, c+from_c_offset, 500
x_range = range(x_start, x_end, x_N)

f(x) = nonlinearity_1(x,a,b,c)
dfdx(x) = gradient(f, x)[1]
f(x::Vector{Float64}) = f(x[1])
dfdx(x::Vector{Float64}) = dfdx(x[1])

y_ = f.(x_range)

#########################################################################################

x_D₁ = Optim.minimizer(Optim.optimize(x->-f(x), [(c+b)/2]))[1]
x_C₂ = Optim.minimizer(Optim.optimize(x-> f(x), [(b+a)/2]))[1]
x_D₃ = Optim.minimizer(Optim.optimize(x->-f(x), [(a+0)/2]))[1]
x_B₁ = -x_D₁
x_A₂ = -x_C₂
x_B₃ = -x_D₃

y_D₁ = f(x_D₁)
y_C₂ = f(x_C₂)
y_D₃ = f(x_D₃)
y_B₁ = f(x_B₁)
y_A₂ = f(x_A₂)
y_B₃ = f(x_B₃)
y_A₁ = y_D₁
y_C₁ = y_B₁
y_D₂ = y_A₂
y_B₂ = y_C₂
y_A₃ = y_D₃
y_C₃ = y_B₃

x_C₁ = Optim.minimizer(Optim.optimize(x->(f(x)-y_C₁)^2, [c], LBFGS()))[1]
x_D₂ = Optim.minimizer(Optim.optimize(x->(f(x)-y_D₂)^2, [b], LBFGS()))[1]
x_C₃ = Optim.minimizer(Optim.optimize(x->(f(x)-y_C₃)^2, [a], LBFGS()))[1]

x_A₁ = -x_C₁
x_B₂ = -x_D₂
x_A₃ = -x_C₃

#########################################################################################

println("A₁ = $((x_A₁, y_A₁))")
println("B₁ = $((x_B₁, y_B₁))")
println("C₁ = $((x_C₁, y_C₁))")
println("D₁ = $((x_D₁, y_D₁))")

println("A₂ = $((x_A₂, y_A₂))")
println("B₂ = $((x_B₂, y_B₂))")
println("C₂ = $((x_C₂, y_C₂))")
println("D₂ = $((x_D₂, y_D₂))")

println("A₃ = $((x_A₃, y_A₃))")
println("B₃ = $((x_B₃, y_B₃))")
println("C₃ = $((x_C₃, y_C₃))")
println("D₃ = $((x_D₃, y_D₃))")


#########################################################################################

fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1], 
    title="", 
    xlabel="u", 
    ylabel="v",
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)
vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

lines!(ax, x_range, y_, color=:blue, label="f(x)")

scatter!(ax, [x_A₁], [y_A₁], label="A₁")
scatter!(ax, [x_B₁], [y_B₁], label="B₁")
scatter!(ax, [x_C₁], [y_C₁], label="C₁")
scatter!(ax, [x_D₁], [y_D₁], label="D₁")

scatter!(ax, [x_A₂], [y_A₂], label="A₂")
scatter!(ax, [x_B₂], [y_B₂], label="B₂")
scatter!(ax, [x_C₂], [y_C₂], label="C₂")
scatter!(ax, [x_D₂], [y_D₂], label="D₂")

scatter!(ax, [x_A₃], [y_A₃], label="A₃")
scatter!(ax, [x_B₃], [y_B₃], label="B₃")
scatter!(ax, [x_C₃], [y_C₃], label="C₃")
scatter!(ax, [x_D₃], [y_D₃], label="D₃")

lines!(ax, [x_B₁, x_C₁], [y_B₁, y_C₁])
lines!(ax, [x_A₁, x_D₁], [y_A₁, y_D₁])
lines!(ax, [x_B₂, x_C₂], [y_B₂, y_C₂])
lines!(ax, [x_A₂, x_D₂], [y_A₂, y_D₂])
lines!(ax, [x_B₃, x_C₃], [y_B₃, y_C₃])
lines!(ax, [x_A₃, x_D₃], [y_A₃, y_D₃])

axislegend(ax, position=:lb) # (l, r, c), (b, t, c)
save(plotsdir("10-nonlinearity_analysis_$(time_ns()).png"), fig, px_per_unit=2)