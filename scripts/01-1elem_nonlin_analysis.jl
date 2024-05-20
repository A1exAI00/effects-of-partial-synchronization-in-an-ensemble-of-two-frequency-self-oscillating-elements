#=
Цель программы: 
Проанализировать нелинейность f(x,a,b,c) при конкретных значениях a,b,c.
Под этим подразумеваются точки и значения экстремумов.

Значения a,b,c: a = 0.32, b = 0.79, c = 1.166
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie
using NLsolve

include(srcdir("article1_module.jl"))
include(srcdir("misc.jl"))
include(srcdir("plotting_tools.jl"))

using .article1

#########################################################################################

x_from_c_offset = 0.06
x_start = - article1.c - x_from_c_offset
x_end = article1.c + x_from_c_offset
x_N = 500
x_range = range(x_start, x_end, x_N)

#########################################################################################

f(x) = article1.f(x)
f(x::Vector{Float64}) = f(x[1])

dfdx(x) = article1.f_deriv(x)
dfdx(x::Vector{Float64}) = dfdx(x[1])

function dfdx!(F, x)
    F[1] = dfdx(x[1])
end

function function_generator(shift)
    function inner_opt!(F, x)
        F[1] = f(x[1]) - shift
    end
    return inner_opt!
end

#########################################################################################

f_res = f.(x_range)

x_D₁ = nlsolve(dfdx!, [(article1.c+article1.b)/2], autodiff=:forward).zero[1]
x_C₂ = nlsolve(dfdx!, [(article1.b+article1.a)/2], autodiff=:forward).zero[1]
x_D₃ = nlsolve(dfdx!, [(article1.a+0)/2], autodiff=:forward).zero[1]
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

x_C₁ = nlsolve(function_generator(y_C₁), [article1.c], autodiff=:forward).zero[1]
x_D₂ = nlsolve(function_generator(y_D₂), [article1.b], autodiff=:forward).zero[1]
x_C₃ = nlsolve(function_generator(y_C₃), [article1.a], autodiff=:forward).zero[1]

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

fig = Figure(size=(700, 500))
ax = beautiful_Axis(fig[1, 1]; title="", xlabel="u", ylabel="v")
vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

lines!(ax, x_range, f_res, color=:cyan, label="f(x)")

lines!(ax, [x_B₁, x_C₁], [y_B₁, y_C₁], color=:blue)
lines!(ax, [x_A₁, x_D₁], [y_A₁, y_D₁], color=:blue)
lines!(ax, [x_B₂, x_C₂], [y_B₂, y_C₂], color=:red)
lines!(ax, [x_A₂, x_D₂], [y_A₂, y_D₂], color=:red)
lines!(ax, [x_B₃, x_C₃], [y_B₃, y_C₃], color=:blue)
lines!(ax, [x_A₃, x_D₃], [y_A₃, y_D₃], color=:blue)

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

limits!(ax, nothing, nothing, -0.2, 0.2)

axislegend(ax, position=:lb) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "01-1elem_nonlin_analysis_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)