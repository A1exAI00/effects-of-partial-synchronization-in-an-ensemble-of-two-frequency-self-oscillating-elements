#=
Цель программы: 
Построить график зависимости средней амплитуды колебаний элементов в цепочке 
в зависимости от параметра d при J=0.

Выбрать НУ соответственно последовательности "1000110", фазы равны 0.
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie
using DifferentialEquations:Rosenbrock23

include(srcdir("article1_module.jl"))
include(srcdir("misc.jl"))
include(srcdir("plotting_tools.jl"))

using .article1

#########################################################################################

N_elements = 7

d_start, d_end, d_N = 0.0, 0.05, 100

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
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

aᵢ_from_d = []
for (i,d) in enumerate(d_range)
    global U₀_tmp
    println(i)
    sol = article1.integrate_multiple_elements(U₀_tmp, t_span, d, N_elements, alg=Rosenbrock23())
    U₀_tmp = sol[:,end]
    uᵢ = [sol[k,:] for k in 1:N_elements]
    vᵢ = [sol[N_elements+k,:] for k in 1:N_elements]

    aᵢ = zeros(N_elements)
    for j in 1:N_elements
        aᵢ[j] = calc_avg_amtlitude(vᵢ[j], uᵢ[j])
    end
    push!(aᵢ_from_d, aᵢ)
end

#########################################################################################

fig = Figure(size=(700, 500))
ax = beautiful_Axis(fig[1, 1], 
    title="Зависимость средних амплитуд элементов цепочки от параметра d\nφ-$φ_mode, t_int=$t_end", 
    xlabel="d", ylabel="⟨aᵢ⟩"
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

for i in 1:N_elements
    lines!(ax, d_range, [aᵢ_from_d[j][i] for j in eachindex(d_range)], label="a$(to_subscript(i))")
end

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "09-Nelem_a(d)_J=0_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)