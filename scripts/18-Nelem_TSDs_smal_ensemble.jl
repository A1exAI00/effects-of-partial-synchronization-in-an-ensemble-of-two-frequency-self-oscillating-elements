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

# function init_cond(i, i_max)
#     if i ≤ i_max/2
#         return 1
#     else
#         if iseven(i)
#             return 0
#         else
#             return 1
#         end
#     end
# end

#########################################################################################

N_elements = 7
d = 0.006
J = 0.0

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end, t_N = 0.0, 3e4, 1000

frames_N = 100
periods_to_avg = 10000

#########################################################################################

init_points = article1.φ_mode_to_init_points(φ_mode, initial_pattern)

u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]

t_spans = [(t_start+t_end*(i-1), t_start+t_end*(i)) for i in 1:frames_N]

U₀_tmp = deepcopy(U₀)

#########################################################################################

Δtₙ = [Float64[] for i in 1:N_elements]

t_timer = time()
Δt_to_avg = Float64[]

for i in eachindex(t_spans)
    global U₀_tmp, t_timer
    Δt_timer_curr = time()-t_timer
    push!(Δt_to_avg, Δt_timer_curr)
    t_timer = time()
    eta_curr = (frames_N - i)*mean(Δt_to_avg)
    println("$i/$(frames_N): Δt=$(round(Δt_timer_curr, digits=5)), eta=$(round(eta_curr, digits=5))")
        

    t_span = t_spans[i]
    
    sol = article1.moded_integrate_multiple_elements(U₀_tmp, t_span, d, N_elements, J; saveat=range(t_span..., t_N))
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
        title="Зависимость средней частоты элементов цепочки от номера элемента; N=$(N_elements), d=$(d), J=$J", 
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
    # limits!(ax1, nothing, nothing, 0.0, 0.0035)
    # limits!(ax2, t_span..., nothing, nothing)

    Colorbar(fig[2, 2], limits = (-π,π), flipaxis = false)

    # axislegend(ax1, position=:rb) # (l, r, c), (b, t, c)
    save_path = plotsdir("21-moded_TSDs_smal_ensemble", "21-moded_TSDs_smal_ensemble_$(lpad(i,3,"0")).png")
    save(save_path, fig, px_per_unit=2)
end

#########################################################################################