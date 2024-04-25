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

v₁⁰, v₂⁰ = 0.0, 0.0
d, J = 0.01, 0.0

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

G_res = [G(u₁, u₂, v₁⁰, v₂⁰, d) for u₁ in u₁_range, u₂ in u₂_range]

#########################################################################################

fig = Figure(size=(700, 700))
ax = beautiful_Axis(fig[1, 1], 
    title="Градиент в системе быстрых движений\nd=$d, v₁⁰=$v₁⁰, v₂⁰=$v₂⁰", 
    xlabel="u₁", ylabel="u₂"
)


contour!(ax, u₁_range, u₂_range, G_res, levels=range(-0.1, maximum(G_res), 100))
vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

limits!(ax, u₁_start, u₁_end, u₂_start, u₂_end)

savepath = plotsdir("tmp", "23-2elem_fast_subsys_PS_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)