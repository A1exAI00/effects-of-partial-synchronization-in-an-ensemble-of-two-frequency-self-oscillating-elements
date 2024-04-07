#=
Цель программы: показать химерное состояние для цепочки из N_elements=1000 элементов.

TODO: уточнить описание программы, что здесь зависимость конечной фазы от параметра d
TODO: пофиксить этот ужас с интексами
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie
using Statistics:mean

include(srcdir("article1_module.jl"))
include(srcdir("misc.jl"))
include(srcdir("plotting_tools.jl"))

using .article1

#########################################################################################

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

N_elements = 100
d = 0.006

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"

t_start, t_end, t_N = 0.0, 1e4, 1000

frames_N = 100

periods_to_avg = 10000

#########################################################################################

initial_pattern = Bool.(init_cond.(1:N_elements, N_elements))
init_points = article1.φ_mode_to_init_points(φ_mode, initial_pattern)

u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]
t_spans = [(t_start+t_end*(i-1), t_start+t_end*(i)) for i in 1:frames_N]

U₀_tmp = deepcopy(U₀)

#########################################################################################

Δtₙ = [Float64[] for i in 1:N_elements]

for i in eachindex(t_spans)
    global U₀_tmp
    println(i)

    t_span = t_spans[i]
    
    sol = article1.integrate_multiple_elements(U₀_tmp, t_span, d, N_elements; saveat=range(t_span..., t_N))
    sol_t = sol.t
    uᵢ = [sol[i,:] for i in 1:N_elements]
    vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

    if i == length(t_spans)
        println("Конечное состояние:\n$(sol[:,end])")
    end

    U₀_tmp = sol[:,end]

    ωᵢ = zeros(N_elements)
    for j in 1:N_elements
        tₙ_curr = passing_through_zero_positive_t(uᵢ[j], sol_t)
        append!(Δtₙ[j], diff(tₙ_curr))
        # ωᵢ[j] = 1/mean(last(Δtₙ[j], periods_to_avg))
        ωᵢ[j] = 1/mean(Δtₙ[j])
    end

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

    ###

    fig = Figure(size=(1000, 700))
    ax1 = beautiful_Axis(fig[1, 1], 
        title="Зависимость средней частоты элементов цепочки от номера элемента; N=$(N_elements), d=$(d), φ-$φ_mode", 
        xlabel="i", ylabel="⟨ωᵢ⟩"
    )
    ax2 = beautiful_Axis(fig[2, 1], 
        title="Зависимость фазы элементов цепочки от времени", 
        xlabel="t", ylabel="i"
    )

    vlines!(ax1, 0.0, color=:black)
    hlines!(ax1, 0.0, color=:black)

    scatter!(ax1, 1:N_elements, ωᵢ, color=:blue)

    heatmap!(ax2, range(t_span..., t_N), 1:N_elements, transpose(φᵢₜ)) 
    limits!(ax1, nothing, nothing, 0.0, 0.0035)
    limits!(ax2, t_span..., nothing, nothing)

    Colorbar(fig[2, 2], limits = (-π,π), flipaxis = false)

    # axislegend(ax1, position=:rb) # (l, r, c), (b, t, c)
    save_path = plotsdir("11-article1_last_plots", "11-article1_last_plot_$(lpad(i,3,"0")).png")
    save(save_path, fig, px_per_unit=2)
end

#########################################################################################