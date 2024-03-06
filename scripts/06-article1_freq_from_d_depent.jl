#=
Цель программы: закономерности смены режимов системв при изменении силы связи d.

Выбрать НУ соответственно последовательности 1000110, фазы распределить случайно.
Для каждого значения d итегрировать систему N_repeat раз, каждый раз записывать 
среднюю амплитуду, частоту и фазу i-го осциллятора в цепочке, потом усреднить.

TODO: переместить функцию calc_avg_freq в misc.jl 
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

#########################################################################################

function calc_avg_freq(d)
    p = (a, b, c, ε, d, N_elements)
    sol = integrate_chain(U₀, p, t_span)
    sol_t = sol.t
    uᵢ = [sol[i,:] for i in 1:N_elements]

    ωᵢ = []
    for i in 1:N_elements
        u = uᵢ[i]
        curr_period = calc_period(u, sol_t)

        push!(ωᵢ, 1/curr_period)
    end
    return ωᵢ
end

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001
N_elements = 7

d_start, d_end, d_N = 0.0, 0.02, 100
d_range = range(d_start, d_end, d_N)

initial_pattern = [true, false, false, false, true, true, false]

u₀ = @. c * initial_pattern + a * !initial_pattern
v₀ = @. 0.0 * ones(Float64, N_elements)

t_start, t_end = 0.0, 5e4


#########################################################################################

U₀ = [u₀..., v₀...]
random_initial_offset = 2 .* (rand(2*N_elements) .- 0.5) .* 0.05
U₀ += random_initial_offset
t_span = [t_start, t_end]

#########################################################################################

ωᵢ_from_d = []
for (i,d) in enumerate(d_range)
    println(i)
    ωᵢ = calc_avg_freq(d)
    push!(ωᵢ_from_d, ωᵢ)
end

#########################################################################################

fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1], 
    title="Зависимость средних частот элементов цепочки от параметра d", 
    xlabel="d", 
    ylabel="ωᵢ", 
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
    lines!(ax, d_range, [ωᵢ_from_d[j][i] for j in eachindex(d_range)], label="ω_$(i)")
end

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("06-article1_freq_from_d_depent_$(time_ns()).png"), fig, px_per_unit=2)
