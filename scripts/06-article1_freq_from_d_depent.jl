#=
Цель программы: закономерности смены режимов системв при изменении силы связи d.

Выбрать НУ соответственно последовательности 1000110, фазы распределить случайно.
Для каждого значения d итегрировать систему N_repeat раз, каждый раз записывать 
среднюю амплитуду, частоту и фазу i-го осциллятора в цепочке, потом усреднить.

TODO: уточнить описание программы, что здесь зависимость средней частоты от параметра d
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

N_elements = 7

d_start, d_end, d_N = 0.0, 0.02, 100

φ_mode = "random" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 1e6

#########################################################################################

d_range = range(d_start, d_end, d_N)

init_points = article1.φ_mode_to_init_points(φ_mode, initial_pattern)

u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]
t_span = [t_start, t_end]

U₀_tmp = deepcopy(U₀)

#########################################################################################

ωᵢ_from_d = []
for (i,d) in enumerate(d_range)
    global U₀_tmp
    println(i)
    sol = article1.integrate_multiple_elements(U₀_tmp, t_span, d, N_elements)
    U₀_tmp = sol[:,end]
    uᵢ = [sol[k,:] for k in 1:N_elements]

    ωᵢ = zeros(N_elements)
    for j in 1:N_elements
        ωᵢ[j] = calc_avg_freq(uᵢ[j], sol.t)
    end
    push!(ωᵢ_from_d, ωᵢ)
end

#########################################################################################

fig = Figure(size=(1000, 700))
ax = beautiful_Axis(fig[1, 1], 
    title="Зависимость средних частот элементов цепочки от параметра d; φ-$φ_mode", 
    xlabel="d", ylabel="⟨ωᵢ⟩"
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

for i in 1:N_elements
    lines!(ax, d_range, [ωᵢ_from_d[j][i] for j in eachindex(d_range)], label="ω_$(i)")
end

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("06-article1_freq_from_d_depent_$(time_ns()).png"), fig, px_per_unit=2)