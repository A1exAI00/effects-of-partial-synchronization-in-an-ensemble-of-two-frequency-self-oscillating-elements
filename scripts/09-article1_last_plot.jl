#=
Цель программы: показать химерное состояние для цепочки из N_elements=1000 элементов.

TODO: Добавить про НУ
TODO: уточнить описание программы, что здесь зависимость конечной фазы от параметра d
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

function box_rand(len)
    return 2 .* (rand(len) .- 0.5)
end

function init_cond(i, i_max)
    if i ≤ i_max/2
        return 1
    else
        if iseven(i)
            return 0
        else
            return 1
        end
    end
end

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001
N_elements = 50
d = 0.02

# d_start, d_end, d_N = 0.0, 0.05, 100
# d_range = range(d_start, d_end, d_N)

initial_pattern = Bool.(init_cond.(1:N_elements, N_elements))
println(initial_pattern)

u₀ = @. c * initial_pattern + a * !initial_pattern
v₀ = @. 0.0 * ones(Float64, N_elements)

t_start, t_end = 0.0, 5e4
t_N = 500

#########################################################################################

U₀ = [u₀..., v₀...] .+ 0.05 .* box_rand(2*N_elements) 
p = (a, b, c, ε, d, N_elements)
t_span = [t_start, t_end]

#########################################################################################

ωᵢ = Float64.(article1_calc_avg_freq(U₀, p, t_span))
# println(size(ωᵢ))
# println(typeof(ωᵢ))
# println(typeof(ωᵢ[1]))
# println(ωᵢ)


sol = integrate_chain(U₀, p, t_span; saveat=range(t_span..., t_N))
sol_t = sol.t
uᵢ = [sol[i,:] for i in 1:N_elements]
vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

φᵢₜ = zeros(N_elements, t_N)
for t in 1:t_N
    curr_φ₁ = 0.0
    for i in 1:N_elements
        if i == 1
            curr_φ₁ = calc_phase(uᵢ[i], sol_t, range(t_span..., t_N)[t])
            φᵢₜ[i, t] = 0.0
        else
            φᵢₜ[i, t] = calc_phase(uᵢ[i], sol_t, range(t_span..., t_N)[t]) - curr_φ₁
        end
    end
end

φᵢₜ = rem2pi.(φᵢₜ, RoundNearest)

#########################################################################################

fig = Figure(size=(1000, 700))
ax1 = Axis(fig[1, 1], 
    title="Зависимость средней частоты элементов цепочки от номера элемента; N=$(N_elements), d=$(d)", 
    xlabel="i", 
    ylabel="⟨ωᵢ⟩", 
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)
ax2 = Axis(fig[2, 1], 
    title="Зависимость фазы элементов цепочки от времени", 
    xlabel="t", 
    ylabel="i", 
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

vlines!(ax1, 0.0, color=:black)
hlines!(ax1, 0.0, color=:black)
# scatter!(ax1, 1:Int(N_elements/2), ωᵢ[1:Int(N_elements/2)], color=:blue, label="Химера")
# scatter!(ax1, Int(N_elements/2+1):N_elements, ωᵢ[Int(N_elements/2+1):N_elements], 
#     color=:red, label="Не химера"
# )
scatter!(ax1, 1:N_elements, ωᵢ, color=:blue, label="Химера")


heatmap!(ax2, range(t_span..., t_N), 1:N_elements, transpose(φᵢₜ)) 

# axislegend(ax1, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("09-article1_last_plot_$(time_ns()).png"), fig, px_per_unit=2)
