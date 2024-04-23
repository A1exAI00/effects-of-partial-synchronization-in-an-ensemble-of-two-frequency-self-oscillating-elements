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

J_start, J_end, J_N = 0.0, 1.05, 1000
d₁, d₂ = 0.05, 0.2

J_range = range(J_start, J_end, J_N)

#########################################################################################

# λ₁_res = article1.λ₁.(J_range)
# λ₂_res = article1.λ₂.(J_range)
λ₃_res₁ = article1.λ₃.(J_range, d₁)
λ₄_res₁ = article1.λ₄.(J_range, d₁)
λ₃_res₂ = article1.λ₃.(J_range, d₂)
λ₄_res₂ = article1.λ₄.(J_range, d₂)

Reλ₃_res₁, Imλ₃_res₁ = real.(λ₃_res₁), imag.(λ₃_res₁)
Reλ₄_res₁, Imλ₄_res₁ = real.(λ₄_res₁), imag.(λ₄_res₁)
Reλ₃_res₂, Imλ₃_res₂ = real.(λ₃_res₂), imag.(λ₃_res₂)
Reλ₄_res₂, Imλ₄_res₂ = real.(λ₄_res₂), imag.(λ₄_res₂)

J_extr = [article1.D₃[1], article1.C₂[1], article1.D₁[1]]


#########################################################################################

fig = Figure(size=(700, 700))
ax1 = beautiful_Axis(fig[1, 1], 
	title="d=$d₁", 
	xlabel="J", ylabel="Re/Im λᵢ"
)
ax2 = beautiful_Axis(fig[2, 1], 
	title="d=$d₂", 
	xlabel="J", ylabel="Re/Im λᵢ"
)

vlines!(ax1, 0.0, color=:black)
hlines!(ax1, 0.0, color=:black)
vlines!.(ax1, J_extr, color=:green, linestyle=:dash)

# lines!(ax1, J_range, Reλ₁_res, label="Reλ₁")
# lines!(ax1, J_range, Reλ₂_res, label="Reλ₂")
# lines!(ax1, J_range, Imλ₁_res, label="Imλ₁")
# lines!(ax1, J_range, Imλ₂_res, label="Imλ₂")

lines!(ax1, J_range, Reλ₃_res₁, label="Reλ₃")
lines!(ax1, J_range, Reλ₄_res₁, label="Reλ₄")
lines!(ax1, J_range, Imλ₃_res₁, label="Imλ₃")
lines!(ax1, J_range, Imλ₄_res₁, label="Imλ₄")

vlines!(ax2, 0.0, color=:black)
hlines!(ax2, 0.0, color=:black)
vlines!.(ax2, J_extr, color=:green, linestyle=:dash)

lines!(ax2, J_range, Reλ₃_res₂, label="Reλ₃")
lines!(ax2, J_range, Reλ₄_res₂, label="Reλ₄")
lines!(ax2, J_range, Imλ₃_res₂, label="Imλ₃")
lines!(ax2, J_range, Imλ₄_res₂, label="Imλ₄")

axislegend(ax1, position=:lt) # (l, r, c), (b, t, c)
axislegend(ax2, position=:lt) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "21-2elem_lambdas_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)