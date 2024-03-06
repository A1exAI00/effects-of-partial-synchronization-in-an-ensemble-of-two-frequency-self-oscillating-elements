#=
Цель программы: изобразить график зависимости реальной и мнимой частей 
характеристических показателей λᵢ, i=1,2 в зависимости от значения параметра J.

Формула для λᵢ была получена аналитически из линеаризованной системы.
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

#########################################################################################

function γ(a,b,c,J)
    a², b², c² = a^2, b^2, c^2
    return a²*b²*c² - 3*(a²*b²+a²*c²+b²*c²)*J^2 + 5*(a²+b²+c²)*J^4 - 7*J^6
end

Reλ₁(γ, ε) = (γ^2-4ε>0) ? 0.5*(γ+sqrt(γ^2-4ε)) : 0.5*γ
Reλ₂(γ, ε) = (γ^2-4ε>0) ? 0.5*(γ-sqrt(γ^2-4ε)) : 0.5*γ
Imλ₁(γ, ε) = (γ^2-4ε>0) ? 0 : 0.5*sqrt(4ε-γ^2)
Imλ₂(γ, ε) = (γ^2-4ε>0) ? 0 : -0.5*sqrt(4ε-γ^2)

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001

J_start, J_end, J_N = 0.0, c-0.07, 1000
J_range = range(J_start, J_end, J_N)

#########################################################################################

γ_ = γ.(a, b, c, J_range)

Reλ₁_ = Reλ₁.(γ_, ε)
Reλ₂_ = Reλ₂.(γ_, ε)
Imλ₁_ = Imλ₁.(γ_, ε)
Imλ₂_ = Imλ₂.(γ_, ε)

f_ = nonlinearity_1.(J_range, a, b, c)

extr_index = find_extremum_index(f_)
J_extr = J_range[extr_index]
f_extr = f_[extr_index]

#########################################################################################

fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1], 
    title="Характеристические показатели в зависимости от параметра J", 
    xlabel="J", 
    ylabel="Re/Im λᵢ",
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

vlines!.(ax, J_extr, color=:green, linestyle=:dash)

lines!(ax, J_range, f_, label="Нелинейность f(U)")

lines!(ax, J_range, Reλ₁_, label="Reλ₁")
lines!(ax, J_range, Reλ₂_, label="Reλ₂")
lines!(ax, J_range, Imλ₁_, label="Imλ₁")
lines!(ax, J_range, Imλ₂_, label="Imλ₂")

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("02-lambdas_$(time_ns()).png"), fig, px_per_unit=2)