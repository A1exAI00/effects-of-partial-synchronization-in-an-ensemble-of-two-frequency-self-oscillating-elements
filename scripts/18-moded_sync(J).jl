#=
Цель скрипта: посмотреть как влияет изменение параметра J на синхронную часть химеры

Для этого выбирается граница по частоте, и те элементы, чья средняя частота 
ниже этой частоты - относятся к синхронной части.
Выводится график зависимости количества элементов в синхронной части к общему 
количеству элементов.
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
d = 0.006
ω_threshold = 0.00235

J_start, J_end, J_N = 0.0, article1.D₃[1], 1000

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 1e5

#########################################################################################

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

n_sync = zeros(length(J_range))

for (i,J) in enumerate(J_range)
    global U₀_tmp, t_timer
    Δt_timer_curr = time()-t_timer
    push!(Δt_to_avg, Δt_timer_curr)
    t_timer = time()
    eta_curr = (J_N - i)*mean(Δt_to_avg)
    println("$i/$(J_N): Δt=$(round(Δt_timer_curr, digits=5)), eta=$(round(eta_curr, digits=5))")
    
    sol = article1.moded_integrate_multiple_elements(U₀_tmp, t_span, d, N_elements, J)
    U₀_tmp = sol[:,end]
    uᵢ = [sol[k,:] for k in 1:N_elements]
    vᵢ = [sol[N_elements+k,:] for k in 1:N_elements]

    ωᵢ = zeros(N_elements)
    for k in 1:N_elements
        ωᵢ[k] = calc_avg_freq(uᵢ[k].-J, sol.t)
    end

    sync_elements_n = sum(ωᵢ .< ω_threshold)
    n_sync[i] = sync_elements_n
end

fig = Figure(size=(1000, 700))
ax = beautiful_Axis(fig[1, 1], 
    title="Зависимость количества элементов в синх. части от параметра J; φ-$φ_mode, d=$d, t_int=$t_end, ω_threshold=$ω_threshold", 
    xlabel="J", ylabel="n"
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

# scatter!(ax, J_range, n_sync, label="")
lines!(ax, J_range, n_sync, label="")

# axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
savepath = plotsdir("18-moded_sync(J)_$(time()).png")
save(savepath, fig, px_per_unit=2)