#=

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

N_elements = 7

d_start, d_end, d_N = 0.0, 0.02, 50

J_start, J_end, J_N = 0.0, article1.D₃[1], 100

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 1e5

#########################################################################################

d_range = range(d_start, d_end, d_N)
J_range = range(J_start, J_end, J_N)

init_points = article1.φ_mode_to_init_points(φ_mode, initial_pattern)
    
u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]
t_span = [t_start, t_end]

U₀_tmp = deepcopy(U₀)

#########################################################################################

t_timer = time()
Δt_to_avg = Float64[]

for (j,d) in enumerate(d_range)
    global U₀_tmp, t_timer
    ωᵢ_from_J = []

    U₀_tmp = deepcopy(U₀)
    for (i,J) in enumerate(J_range)

        Δt_timer_curr = time()-t_timer
        push!(Δt_to_avg, Δt_timer_curr)
        t_timer = time()
        eta_curr = (d_N*J_N - (j*J_N + i))*mean(Δt_to_avg)
        println("$j/$(d_N), $i/$(J_N): Δt=$(round(Δt_timer_curr, digits=5)), eta=$(round(eta_curr, digits=5))")
        
        sol = article1.moded_integrate_multiple_elements(U₀_tmp, t_span, d, N_elements, J)
        U₀_tmp = sol[:,end]
        uᵢ = [sol[k,:] for k in 1:N_elements]
        vᵢ = [sol[N_elements+k,:] for k in 1:N_elements]

        ωᵢ = zeros(N_elements)
        for k in 1:N_elements
            ωᵢ[k] = calc_avg_freq(uᵢ[k].-J, sol.t)
        end

        push!(ωᵢ_from_J, ωᵢ)
    end

    #########################################################################################

    fig = Figure(size=(1000, 700))
    ax = beautiful_Axis(fig[1, 1], 
        title="Зависимость средней частоты элементов цепочки от параметра J; φ-$φ_mode, d=$d, t_int=$t_end", 
        xlabel="J", ylabel="ωᵢ"
    )

    vlines!(ax, 0.0, color=:black)
    hlines!(ax, 0.0, color=:black)

    for i in 1:N_elements
        scatter!(ax, J_range, [ωᵢ_from_J[k][i] for k in eachindex(J_range)], label="φ_$(i) - φ₅")
    end

    axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
    savepath = plotsdir("15-moded_freq_from_J", "15-moded_freq_from_J_$(lpad(j,3,"0")).png")
    save(savepath, fig, px_per_unit=2)
end