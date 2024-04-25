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

v₁⁰_11, v₂⁰_11, d_11 = 0.0, 0.0, 0.00
v₁⁰_12, v₂⁰_12, d_12 = 0.0, 0.02, 0.00
v₁⁰_21, v₂⁰_21, d_21 = 0.0, 0.0, 0.02
v₁⁰_22, v₂⁰_22, d_22 = 0.0, 0.02, 0.02
v₁⁰_31, v₂⁰_31, d_31 = 0.0, 0.0, 0.03
v₁⁰_32, v₂⁰_32, d_32 = 0.0, 0.02, 0.03

offset = 0.2
u₁_start, u₁_end, u₁_N = -article1.c-offset, article1.c+offset, 1000
u₂_start, u₂_end, u₂_N = -article1.c-offset, article1.c+offset, 1000

#########################################################################################

u₁_range = range(u₁_start, u₁_end, u₁_N)
u₂_range = range(u₂_start, u₂_end, u₂_N)

function G(u₁, u₂, v₁⁰, v₂⁰, d)
    a,b,c = article1.a, article1.b, article1.c
    a²,b²,c² = a^2, b^2, c^2
    int_f(x) = 1/8*x^8 - 1/6*(a²+b²+c²)*x^6 + 1/4*(a²*b²+a²*c²+b²*c²)*x^4 - 1/2*(a²*b²*c²)*x^2

    G_curr = - int_f(u₁) - int_f(u₂) + v₁⁰*u₁ + v₂⁰*u₂ + d/2*(u₁-u₂)^2
    return G_curr
end

G_res_11 = [G(u₁, u₂, v₁⁰_11, v₂⁰_11, d_11) for u₁ in u₁_range, u₂ in u₂_range]
G_res_12 = [G(u₁, u₂, v₁⁰_12, v₂⁰_12, d_12) for u₁ in u₁_range, u₂ in u₂_range]
G_res_21 = [G(u₁, u₂, v₁⁰_21, v₂⁰_21, d_21) for u₁ in u₁_range, u₂ in u₂_range]
G_res_22 = [G(u₁, u₂, v₁⁰_22, v₂⁰_22, d_22) for u₁ in u₁_range, u₂ in u₂_range]
G_res_31 = [G(u₁, u₂, v₁⁰_31, v₂⁰_31, d_31) for u₁ in u₁_range, u₂ in u₂_range]
G_res_32 = [G(u₁, u₂, v₁⁰_32, v₂⁰_32, d_32) for u₁ in u₁_range, u₂ in u₂_range]

#########################################################################################

fig = Figure(size=(700, 1000))
ax_11 = beautiful_Axis(fig[1, 1], title="d=$d_11, v₁⁰=$v₁⁰_11, v₂⁰=$v₂⁰_11", xlabel="u₁", ylabel="u₂")
ax_12 = beautiful_Axis(fig[1, 2], title="d=$d_12, v₁⁰=$v₁⁰_12, v₂⁰=$v₂⁰_12", xlabel="u₁", ylabel="u₂")
ax_21 = beautiful_Axis(fig[2, 1], title="d=$d_21, v₁⁰=$v₁⁰_21, v₂⁰=$v₂⁰_21", xlabel="u₁", ylabel="u₂")
ax_22 = beautiful_Axis(fig[2, 2], title="d=$d_22, v₁⁰=$v₁⁰_22, v₂⁰=$v₂⁰_22", xlabel="u₁", ylabel="u₂")
ax_31 = beautiful_Axis(fig[3, 1], title="d=$d_31, v₁⁰=$v₁⁰_31, v₂⁰=$v₂⁰_31", xlabel="u₁", ylabel="u₂")
ax_32 = beautiful_Axis(fig[3, 2], title="d=$d_32, v₁⁰=$v₁⁰_32, v₂⁰=$v₂⁰_32", xlabel="u₁", ylabel="u₂")

contour!(ax_11, u₁_range, u₂_range, G_res_11, levels=range(-0.02, maximum(G_res_11), 100))
contour!(ax_12, u₁_range, u₂_range, G_res_12, levels=range(-0.02, maximum(G_res_12), 100))
contour!(ax_21, u₁_range, u₂_range, G_res_21, levels=range(-0.02, maximum(G_res_21), 100))
contour!(ax_22, u₁_range, u₂_range, G_res_22, levels=range(-0.02, maximum(G_res_22), 100))
contour!(ax_31, u₁_range, u₂_range, G_res_31, levels=range(-0.02, maximum(G_res_31), 100))
contour!(ax_32, u₁_range, u₂_range, G_res_32, levels=range(-0.02, maximum(G_res_32), 100))
vlines!.([ax_11, ax_12, ax_21, ax_22, ax_31, ax_32], 0.0, color=:black)
hlines!.([ax_11, ax_12, ax_21, ax_22, ax_31, ax_32], 0.0, color=:black)

limits!.([ax_11, ax_12, ax_21, ax_22, ax_31, ax_32], u₁_start, u₁_end, u₂_start, u₂_end)

savepath = plotsdir("tmp", "24-2elem_fastsys_PS_LaTeX_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)