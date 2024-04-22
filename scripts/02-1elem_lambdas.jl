#=
Цель программы: изобразить график зависимости реальной и мнимой частей 
характеристических показателей λᵢ, i=1,2 в зависимости от значения параметра J.

Формула для λᵢ была получена аналитически из линеаризованной системы.
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie

include(srcdir("article1_module.jl"))
include(srcdir("misc.jl"))
include(srcdir("plotting_tools.jl"))

using .article1

#########################################################################################

J_start, J_end, J_N = 0.0, 1.05, 1000
J_range = range(J_start, J_end, J_N)

#########################################################################################

Reλ₁_res = article1.Reλ₁.(J_range)
Reλ₂_res = article1.Reλ₂.(J_range)
Imλ₁_res = article1.Imλ₁.(J_range)
Imλ₂_res = article1.Imλ₂.(J_range)

J_extr = [article1.D₃[1], article1.C₂[1], article1.D₁[1]]

L₁ = article1.first_Lyapunov_quantity.(J_range)

#########################################################################################

fig = Figure(size=(700, 700))
ax1 = beautiful_Axis(fig[1, 1], 
	title="Зависимость характеристических показателей λᵢ от параметра J", 
	xlabel="J", ylabel="Re/Im λᵢ"
)
ax2 = beautiful_Axis(fig[2, 1], 
	title="Зависимость первой ляпуновской величины L₁ от параметра J", 
	xlabel="J", ylabel="L₁"
)

vlines!(ax1, 0.0, color=:black)
hlines!(ax1, 0.0, color=:black)

vlines!.(ax1, J_extr, color=:green, linestyle=:dash)

lines!(ax1, J_range, Reλ₁_res, label="Reλ₁")
lines!(ax1, J_range, Reλ₂_res, label="Reλ₂")
lines!(ax1, J_range, Imλ₁_res, label="Imλ₁")
lines!(ax1, J_range, Imλ₂_res, label="Imλ₂")


vlines!(ax2, 0.0, color=:black)
hlines!(ax2, 0.0, color=:black)

vlines!.(ax2, J_extr, color=:green, linestyle=:dash)

this_label = "L₁(J=J₁) = $(article1.first_Lyapunov_quantity(J_extr[1]))\n
L₁(J=J₂) = $(article1.first_Lyapunov_quantity(J_extr[2]))\n
L₁(J=J₃) = $(article1.first_Lyapunov_quantity(J_extr[3]))"
lines!(ax2, J_range, L₁, label=this_label)

axislegend(ax1, position=:lt) # (l, r, c), (b, t, c)
axislegend(ax2, position=:lt) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "02-1elem_lambdas_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)