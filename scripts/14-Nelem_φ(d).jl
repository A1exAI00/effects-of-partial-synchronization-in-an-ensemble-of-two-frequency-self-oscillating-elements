#=
Цель программы: 
Построить график зависимости мгновенной фазы колебаний элементов относительно фазы 
пятого элемента в цепочке в зависимости от параметра d при J≠0.
Для каждого значения параметра J строится отдельный график.

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

savedir = joinpath("tmp", "14-Nelem_φ(d)")
if !(ispath(savedir)) mkpath(plotsdir(savedir)) end

#########################################################################################

N_elements = 7

d_start, d_end, d_N = 0.0, 0.05, 100

J_start, J_end, J_N = 0.0, article1.D₃[1], 50

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 1e4

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

for j in eachindex(J_range)
    global U₀_tmp, t_timer
    J = J_range[j]
    φᵢ_from_d = []

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

        φᵢ = zeros(N_elements)
        for j in 1:N_elements
            φᵢ[j] = calc_phase(uᵢ[j].-J, sol.t, sol.t[end])
        end

        φ₅ = φᵢ[5]
        φᵢ = rem2pi.(φᵢ .- φ₅, RoundNearest)
        push!(φᵢ_from_d, φᵢ)
    end

    #########################################################################################

    fig = Figure(size=(1000, 700))
    ax = beautiful_Axis(fig[1, 1], 
        title="Зависимость конечной фазы элементов цепочки от параметра d; φ-$φ_mode, J=$J, t_int=$t_end", 
        xlabel="d", ylabel="φᵢ"
    )

    vlines!(ax, 0.0, color=:black)
    hlines!(ax, 0.0, color=:black)

    for i in 1:N_elements
        scatter!(ax, d_range, [φᵢ_from_d[j][i] for j in eachindex(d_range)], label="φ_$(i) - φ₅")
    end

    limits!(ax, nothing, nothing, -π, π)

    axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
    
    savepath = plotsdir(savedir, "14-Nelem_φ(d)_$(lpad(j,3,"0")).png")
    save(savepath, fig, px_per_unit=2)
end