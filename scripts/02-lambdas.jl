#=
Цель программы: изобразить график зависимости реальной и мнимой частей 
характеристических показателей λᵢ, i=1,2 в зависимости от значения параметра J.

Формула для λᵢ была получена аналитически из линеаризованной системы.

TODO: заменить нахождение экстремума нелинейности f(u) за более точный метод
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001

J_start, J_end, J_N = 0.0, 1.05, 1000000
J_range = range(J_start, J_end, J_N)

#########################################################################################

γ_ = ∂F₁_∂U_at_eq.(J_range)

Reλ₁_ = Reλ₁.(γ_, ε)
Reλ₂_ = Reλ₂.(γ_, ε)
Imλ₁_ = Imλ₁.(γ_, ε)
Imλ₂_ = Imλ₂.(γ_, ε)

f_ = nonlinearity_1.(J_range, a, b, c)

extr_index = find_extremum_index(f_)
J_extr = J_range[extr_index]
f_extr = f_[extr_index]

L₁ = first_Lyapunov_quantity.(J_range)

#########################################################################################

fig = Figure(size=(1000, 700))
ax1 = Axis(fig[1, 1], 
    title="Зависимость характеристических показателей λᵢ от параметра J", 
    xlabel="J", 
    ylabel="Re/Im λᵢ",
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)
ax2 = Axis(fig[2, 1], 
    title="Зависимость первой ляпуновской величины L₁ от параметра J", 
    xlabel="J", 
    ylabel="L₁",
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

vlines!(ax1, 0.0, color=:black)
hlines!(ax1, 0.0, color=:black)

vlines!.(ax1, J_extr, color=:green, linestyle=:dash)

# lines!(ax1, J_range, f_, label="Нелинейность f(U)")
lines!(ax1, J_range, Reλ₁_, label="Reλ₁")
lines!(ax1, J_range, Reλ₂_, label="Reλ₂")
lines!(ax1, J_range, Imλ₁_, label="Imλ₁")
lines!(ax1, J_range, Imλ₂_, label="Imλ₂")


vlines!(ax2, 0.0, color=:black)
hlines!(ax2, 0.0, color=:black)

vlines!.(ax2, J_extr, color=:green, linestyle=:dash)
lines!(ax2, J_range, L₁, label="L₁(J=J₁) = $(first_Lyapunov_quantity(J_extr[1]))\nL₁(J=J₂) = $(first_Lyapunov_quantity(J_extr[2]))\nL₁(J=J₃) = $(first_Lyapunov_quantity(J_extr[3]))")

axislegend(ax1, position=:lt) # (l, r, c), (b, t, c)
axislegend(ax2, position=:lt) # (l, r, c), (b, t, c)
save(plotsdir("02-lambdas_$(time_ns()).png"), fig, px_per_unit=2)