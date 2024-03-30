#=
Цель скрипта: построить несколько зависимостей конечной фазы от величины 
параметра J для нескольких параметров d.
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

d_start, d_end, d_N = 0.0, 0.05, 50

J_start, J_end, J_N = 0.0, article1.a, 100

φ_mode = "random" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 1e5

#########################################################################################

d_range = range(d_start, d_end, d_N)
J_range = range(J_start, J_end, J_N)

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

U₀_tmp = deepcopy(U₀)

#########################################################################################

for (j,d) in enumerate(d_range)
    global U₀_tmp
    φᵢ_from_J = []

    U₀_tmp = deepcopy(U₀)
    for (i,J) in enumerate(J_range)
        println(i)
        
        sol = article1.moded_integrate_multiple_elements(U₀_tmp, t_span, d, N_elements, J)
        U₀_tmp = sol[:,end]
        uᵢ = [sol[k,:] for k in 1:N_elements]
        vᵢ = [sol[N_elements+k,:] for k in 1:N_elements]

        φᵢ = zeros(N_elements)
        for j in 1:N_elements
            φᵢ[j] = calc_phase(uᵢ[j].-J, sol.t, sol.t[end])
        end

        φ₅ = φᵢ[5]
        φᵢ = rem2pi.(φᵢ .- φ₅, RoundNearest)
        push!(φᵢ_from_J, φᵢ)
    end

    #########################################################################################

    fig = Figure(size=(1000, 700))
    ax = beautiful_Axis(fig[1, 1], 
        title="Зависимость конечной фазы элементов цепочки от параметра J; φ-$φ_mode, d=$d", 
        xlabel="J", ylabel="φᵢ"
    )

    vlines!(ax, 0.0, color=:black)
    hlines!(ax, 0.0, color=:black)

    for i in 1:N_elements
        scatter!(ax, J_range, [φᵢ_from_J[j][i] for j in eachindex(J_range)], label="φ_$(i) - φ₅")
    end

    limits!(ax, nothing, nothing, -π, π)

    axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
    savepath = plotsdir("13-moded_phase_from_J", "13-moded_phase_from_J_$(lpad(j,3,"0")).png")
    save(savepath, fig, px_per_unit=2)
end