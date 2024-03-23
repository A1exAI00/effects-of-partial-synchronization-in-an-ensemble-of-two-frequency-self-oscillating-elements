#=
Цель программы: построить эпюры каждого элемента цепочки.

Зафиксировать все параметры.
Выбрать НУ соответственно последовательности 1000110, фазы распределить случайно.
Построить график uᵢ(t).
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

d, N_elements = 0.02, 7

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 5e3

#########################################################################################

if φ_mode == "random"
    init_points = article1.initial_random_phase.(initial_pattern)
elseif φ_mode == "zero"
    init_points = article1.initial_zero_phase.(initial_pattern)
elseif φ_mode == "синфазно"
    init_points = [article1.initial_ring_phase(initial_pattern[i], i) for i in eachindex(initial_pattern)]
elseif φ_mode == "противофазно"
    init_points = [article1.initial_asinphase(initial_pattern[i], i) for i in eachindex(initial_pattern)]
else
    error("Неверно выбран параметр φ_mode=$φ_mode")
end

u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]
t_span = [t_start, t_end]

#########################################################################################

sol = article1.integrate_multiple_elements(U₀, t_span, d, N_elements)
sol_t = sol.t
uᵢ = [sol[i,:] for i in 1:N_elements]
vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

#########################################################################################

fig = Figure(size=(1000, 700))
ax = beautiful_Axis(fig[1, 1], 
	title="Временные реализации для цепочки элементов; d=$d, φ-$φ_mode", 
	xlabel="t", ylabel="uᵢ"
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

for i in 1:N_elements
    lines!(ax, sol_t, uᵢ[i], label="№$(i)")
end

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("05-article1_chain_$(time_ns()).png"), fig, px_per_unit=2)