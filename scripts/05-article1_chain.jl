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

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001
d = 0.02
N_elements = 7

N_repeat = 10

initial_pattern = [true, false, false, false, true, true, false]

u₀ = @. c * initial_pattern + a * !initial_pattern
v₀ = @. 0.0 * ones(Float64, N_elements)

t_start, t_end = 0.0, 5e3

#########################################################################################

U₀ = [u₀..., v₀...]
random_initial_offset = 2 .* (rand(2*N_elements) .- 0.5) .* 0.1
U₀ += random_initial_offset
p = (a, b, c, ε, d, N_elements)
t_span = [t_start, t_end]

#########################################################################################

sol = integrate_chain(U₀, p, t_span)
sol_t = sol.t
uᵢ = [sol[i,:] for i in 1:N_elements]
vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

#########################################################################################

fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1], 
    title="Временные реализации для цепочки элементов", 
    xlabel="t", 
    ylabel="uᵢ",
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

for i in 1:N_elements
    lines!(ax, sol_t, uᵢ[i], label="№$(i)")
end

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("05-article1_chain_$(time_ns()).png"), fig, px_per_unit=2)