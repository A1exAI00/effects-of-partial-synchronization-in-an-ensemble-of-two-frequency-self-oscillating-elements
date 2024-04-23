#=
Цель программы: 
Изобразить график зависимости реальной и мнимой частей хар. показателей 
λᵢ, i=(1..4) в зависимости от значения параметра J для системы из двух элементов.

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

savedir = joinpath("tmp", "22-2elem_lambdas_multiple")
if !(ispath(savedir)) mkpath(plotsdir(savedir)) end

#########################################################################################

J_start, J_end, J_N = 0.0, 1.05, 1000
d_start, d_end, d_N = 0.0, 1.05, 101

#########################################################################################

J_range = range(J_start, J_end, J_N)
d_range = range(d_start, d_end, d_N)

J_extr = [article1.D₃[1], article1.C₂[1], article1.D₁[1]]

#########################################################################################

function gen_λ(d)
	λ₁_res = article1.λ₁.(J_range)
	λ₂_res = article1.λ₂.(J_range)
	λ₃_res = article1.λ₃.(J_range, d)
	λ₄_res = article1.λ₄.(J_range, d)

	Reλ₁_res, Imλ₁_res = real.(λ₁_res), imag.(λ₁_res)
	Reλ₂_res, Imλ₂_res = real.(λ₂_res), imag.(λ₂_res)
	Reλ₃_res, Imλ₃_res = real.(λ₃_res), imag.(λ₃_res)
	Reλ₄_res, Imλ₄_res = real.(λ₄_res), imag.(λ₄_res)

	return Reλ₁_res, Imλ₁_res, Reλ₂_res, Imλ₂_res, Reλ₃_res, Imλ₃_res, Reλ₄_res, Imλ₄_res
end

#########################################################################################

for (i,d) in enumerate(d_range)

	Reλ₁_res, Imλ₁_res, Reλ₂_res, Imλ₂_res, Reλ₃_res, Imλ₃_res, Reλ₄_res, Imλ₄_res = gen_λ(d)

	fig = Figure(size=(700, 700))
	ax1 = beautiful_Axis(fig[1, 1], 
		title="Зависимость характеристических показателей λᵢ от параметра J\nd=$d", 
		xlabel="J", ylabel="Re/Im λᵢ"
	)
	ax2 = beautiful_Axis(fig[2, 1], 
		title="", 
		xlabel="J", ylabel="Re/Im λᵢ"
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

	lines!(ax2, J_range, Reλ₃_res, label="Reλ₃")
	lines!(ax2, J_range, Reλ₄_res, label="Reλ₄")
	lines!(ax2, J_range, Imλ₃_res, label="Imλ₃")
	lines!(ax2, J_range, Imλ₄_res, label="Imλ₄")

	axislegend(ax1, position=:lt) # (l, r, c), (b, t, c)
	axislegend(ax2, position=:lt) # (l, r, c), (b, t, c)

	savepath = plotsdir(savedir, "22-2elem_lambdas_multiple_$(lpad(i, 3, "0")).png")
	save(savepath, fig, px_per_unit=2)
end

