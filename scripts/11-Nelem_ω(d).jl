#=
Цель программы: 
Построить график зависимости средней частоты колебаний элементов в цепочке 
в зависимости от параметра d при J≠0.

Выбрать НУ соответственно последовательности "1000110", фазы равны 0.
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

savedir = joinpath("tmp", "11-Nelem_ω(d)")
if !(ispath(savedir)) mkpath(plotsdir(savedir)) end

#########################################################################################

N_elements = 7

d_start, d_end, d_N = 0.0, 0.02, 100

J_start, J_end, J_N = 0.0, article1.D₃[1], 50

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

for (j,J) in enumerate(J_range)
    global U₀_tmp, t_timer
    ωᵢ_from_d = []

    U₀_tmp = deepcopy(U₀)
    for (i,d) in enumerate(d_range)

        Δt_timer_curr = time()-t_timer
        push!(Δt_to_avg, Δt_timer_curr)
        t_timer = time()
        eta_curr = (d_N*J_N - (j*d_N + i))*mean(Δt_to_avg)
        println("$j/$(J_N), $i/$(d_N): Δt=$(round(Δt_timer_curr, digits=5)), eta=$(round(eta_curr, digits=5))")
        
        sol = article1.moded_integrate_multiple_elements(U₀_tmp, t_span, d, N_elements, J)
        U₀_tmp = sol[:,end]
        uᵢ = [sol[k,:] for k in 1:N_elements]
        vᵢ = [sol[N_elements+k,:] for k in 1:N_elements]

        ωᵢ = zeros(N_elements)
        for k in 1:N_elements
            ωᵢ[k] = calc_avg_freq(uᵢ[k].-J, sol.t)
        end

        push!(ωᵢ_from_d, ωᵢ)
    end

    #########################################################################################

    fig = Figure(size=(700, 500))
    ax = beautiful_Axis(fig[1, 1], 
        title="Зависимость средней частоты элементов цепочки от параметра d\nφ-$φ_mode, J=$J, t_int=$t_end", 
        xlabel="d", ylabel="⟨ωᵢ⟩"
    )

    vlines!(ax, 0.0, color=:black)
    # hlines!(ax, 0.0, color=:black)

    for i in 1:N_elements
        label_curr = "⟨ω$(to_subscript(i))⟩"
        lines!(ax, d_range, [ωᵢ_from_d[k][i] for k in eachindex(d_range)], label=label_curr)
    end

    axislegend(ax, position=:rb) # (l, r, c), (b, t, c)

    savepath = plotsdir(savedir, "11-Nelem_ω(d)_$(lpad(j,3,"0")).png")
    save(savepath, fig, px_per_unit=2)
end