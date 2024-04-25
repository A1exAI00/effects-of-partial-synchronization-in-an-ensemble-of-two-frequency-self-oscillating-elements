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

d = 0.006
J = 0.17

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start_preintegrate, t_end_preintegrate = 0.0, 1e6
t_start_PHS, t_end_PHS, t_N_PHS = t_end_preintegrate, t_end_preintegrate+4000, 10000

#########################################################################################

init_points = article1.φ_mode_to_init_points(φ_mode, initial_pattern)
    
u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]

t_span_preintegrate = [t_start_preintegrate, t_end_preintegrate]
t_span_PHS = [t_start_PHS, t_end_PHS]

U₀_tmp = deepcopy(U₀)

#########################################################################################

sol_J = article1.moded_integrate_multiple_elements(U₀, t_span_preintegrate, d, N_elements, J)
new_U₀ = sol_J[:,end]

sol_J = article1.moded_integrate_multiple_elements(new_U₀, t_span_PHS, d, N_elements, J, saveat=range(t_span_PHS..., t_N_PHS))
uᵢ_J = [sol_J[k,:] for k in 1:N_elements]
vᵢ_J = [sol_J[N_elements+k,:] for k in 1:N_elements]

ωᵢ = zeros(N_elements)
for k in 1:N_elements
    ωᵢ[k] = calc_avg_freq(sol_J[k,:].-J, sol_J.t)
end

φᵢₜ_J = zeros(N_elements, t_N_PHS)
for t in 1:t_N_PHS
    for i in 1:N_elements
        φᵢₜ_J[i, t] = calc_phase(uᵢ_J[i], sol_J.t, range(t_span_PHS..., t_N_PHS)[t])
    end
end
φᵢₜ_J = rem2pi.(φᵢₜ_J, RoundNearest)

#########################################################################################

fig = Figure(size=(700, 700))

ax_upper = beautiful_Axis(fig[1,1:2],
    title="РАспределение средних частот по номеру элемента\nd=$d, J=$J, ", 
    xlabel="i", ylabel="⟨ωᵢ⟩"
)

ax_11 = beautiful_Axis(fig[2,1],
    title="Парциальный фазов. портрет\nd=$d, J=$J", xlabel="uᵢ", ylabel="vᵢ"
)
ax_12 = beautiful_Axis(fig[2,2],
    title="Простр.-времен. диагр\nd=$d, J=$J", xlabel="t", ylabel="i"
)

scatter!(ax_upper, 1:N_elements, ωᵢ)

vlines!(ax_11, 0.0, color=:black)
hlines!(ax_11, 0.0, color=:black)

for i in 1:N_elements
    tail = round(Int,t_N_PHS/2)
    lines!(ax_11, uᵢ_J[i][end-tail:end], vᵢ_J[i][end-tail:end])
    # lines!(ax_11, uᵢ_J[i], vᵢ_J[i])
end

heatmap!(ax_12, range(t_span_PHS..., t_N_PHS), 1:N_elements, transpose(φᵢₜ_J))
limits!(ax_12, t_span_PHS[2]-4000, t_span_PHS[2], nothing, nothing)

# axislegend(ax, position=:rb) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "20-Nelem_PHS_TSD_LaTeX_$(time_ns()).png")
save(savepath, fig, px_per_unit=3)