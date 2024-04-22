#=
Цель скрипта: посмотреть как будут выглядеть фазовые траектории цепочки элементов
на фазовом пространстве при изменении параметра J.
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

# x range for nonlinearity
x_from_c_offset = 0.05
Δx = 1e-5

N_elements = 7
d = 0.01
J = 0.4

φ_mode = "zero" # "random", "zero", "синфазно", "противофазно"
initial_pattern = [true, false, false, false, true, true, false]

t_start, t_end = 0.0, 1e3

#########################################################################################

x_start = -article1.c - x_from_c_offset
x_end = article1.c - x_from_c_offset
x_N = round(Int,(x_end - x_start)/Δx)
x_range = range(x_start, x_end, x_N)

# nonlinearity
f_res = article1.f.(x_range)
x_extr = [article1.D₃[1], article1.C₂[1], article1.D₁[1]]
y_extr = [article1.D₃[2], article1.C₂[2], article1.D₁[2]]

init_points = article1.φ_mode_to_init_points(φ_mode, initial_pattern)
u₀ = [init_points[i][1] for i in eachindex(init_points)]
v₀ = [init_points[i][2] for i in eachindex(init_points)]
U₀ = [u₀..., v₀...]

t_span = (t_start, t_end)

#########################################################################################

sol = article1.moded_integrate_multiple_elements(U₀, t_span, d, N_elements, J)
uᵢ = [sol[k,:] for k in 1:N_elements]
vᵢ = [sol[N_elements+k,:] for k in 1:N_elements]

#########################################################################################

fig = Figure(size=(700, 500))
ax = beautiful_Axis(fig[1, 1], 
    title="Парциальные фазовые траектории элементов цепочки\nφ_mode=$φ_mode, d=$d, J=$J", 
    xlabel="uᵢ", ylabel="vᵢ"
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

for i in 1:N_elements
    lines!(ax, uᵢ[i], vᵢ[i], label="№$(i)")
end

axislegend(ax, position=:rt) # (l, r, c), (b, t, c)
save(plotsdir("16-moded_phase_spase_$(time_ns()).png"), fig, px_per_unit=2)