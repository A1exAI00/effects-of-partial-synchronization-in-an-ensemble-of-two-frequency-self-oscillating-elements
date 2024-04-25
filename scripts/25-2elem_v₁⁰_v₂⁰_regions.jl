#=

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

v_crit_1, v_crit_2, v_crit_3 = article1.D₃[2], article1.D₂[2], article1.D₁[2]

#########################################################################################

fig = Figure(size=(700, 700))
ax = beautiful_Axis(fig[1, 1], 
    title="d=0", 
    xlabel="v₁⁰", ylabel="v₂⁰"
)

vlines!.(ax, 0.0, color=:black)
hlines!.(ax, 0.0, color=:black)

vlines!.(ax, [v_crit_1, v_crit_2, v_crit_3, -v_crit_1, -v_crit_2, -v_crit_3], color=:green, linestyle=:dash)
hlines!.(ax, [v_crit_1, v_crit_2, v_crit_3, -v_crit_1, -v_crit_2, -v_crit_3], color=:green, linestyle=:dash)

offset = 0.01
limits!.(ax, -v_crit_3-offset, v_crit_3+offset, -v_crit_3-offset, v_crit_3+offset)

savepath = plotsdir("tmp", "25-2elem_v₁⁰_v₂⁰_regions_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)