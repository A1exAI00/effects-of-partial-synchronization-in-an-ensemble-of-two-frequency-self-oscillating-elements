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

u₁₀, v₁₀ = article1.D₃
u₂₀, v₂₀ = article1.A₁

d, J = 0.01, 0.0

t_start, t_end = 0.0, 10000.0

offset = 0.2
u₁_start, u₁_end, u₁_N = -article1.c-offset, article1.c+offset, 1000
u₂_start, u₂_end, u₂_N = -article1.c-offset, article1.c+offset, 1000

#########################################################################################

U₀ = [u₁₀, u₂₀, v₁₀, v₂₀]

t_span = [t_start, t_end]

u₁_range = range(u₁_start, u₁_end, u₁_N)
u₂_range = range(u₂_start, u₂_end, u₂_N)

#########################################################################################

function is_in_fast_movement_lower(x, y)
    if (article1.A₃[1]<x<article1.D₃[1]) && (y>article1.D₃[2])
        return true
    elseif (article1.B₃[1]<x<article1.C₃[1]) && (y<article1.C₃[2])
        return true
    end
    return false
end

function is_in_fast_movement_upper(x, y)
    if (article1.A₁[1]<x<article1.D₁[1]) && (y>article1.D₁[2])
        return true
    elseif (article1.B₁[1]<x<article1.C₁[1]) && (y<article1.C₁[2])
        return true
    end
    return false
end

function is_in_fast_movement(x, y)
    return is_in_fast_movement_lower(x, y) || is_in_fast_movement_upper(x, y)
end

function filter_slow_movement(sol_x1, sol_x2, sol_y1, sol_y2)
    sol_x1_slow, sol_y1_slow = Float64[], Float64[]
    sol_x2_slow, sol_y2_slow = Float64[], Float64[]

    for i in eachindex(sol_x1)
        x1, y1 = sol_x1[i], sol_y1[i]
        x2, y2 = sol_x2[i], sol_y2[i]
        if !(is_in_fast_movement(x1, y1)) && !(is_in_fast_movement(x2, y2))
            push!(sol_x1_slow, x1)
            push!(sol_y1_slow, y1)
            push!(sol_x2_slow, x2)
            push!(sol_y2_slow, y2)
        end
    end

    return sol_x1_slow, sol_x2_slow, sol_y1_slow, sol_y2_slow
end

function Δ_det(u₁, u₂, d)
    f_deriv = article1.f_deriv
    return f_deriv(u₁)*f_deriv(u₂) - d*(f_deriv(u₁)+f_deriv(u₂))
end

function Sₐ(u₁, u₂, d)
    if (Δ_det(u₁, u₂, d) > 0)
        return 1.0
    end
    return 0.0
end

#########################################################################################

sol_preintegrate = article1.moded_integrate_two_elements(U₀, t_span, d, J)
U₀ = sol_preintegrate[:,end]

sol = article1.moded_integrate_two_elements(U₀, t_span, d, J)

sol_u₁ = sol[1,:]
sol_u₂ = sol[2,:]
sol_v₁ = sol[3,:]
sol_v₂ = sol[4,:]

sol_u₁_slow, sol_u₂_slow, sol_v₁_slow, sol_v₂_slow = filter_slow_movement(sol_u₁, sol_u₂, sol_v₁, sol_v₂)

Sₐ_res = [Sₐ(u₁, u₂, d) for u₁ in u₁_range, u₂ in u₂_range]

#########################################################################################

fig = Figure(size=(700, 700))
ax = beautiful_Axis(fig[1, 1], 
    title="Фазовый портрет подсистемы медленных движений\nd=$d, J=$J", 
    xlabel="u₁", ylabel="u₂"
)

contour!(ax, u₁_range, u₂_range, Sₐ_res, levels=10)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

# lines!(ax, sol_u₁_slow, sol_u₂_slow, alpha=0.1)
scatter!(ax, sol_u₁_slow, sol_u₂_slow)

# axislegend(ax, position=:rb) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "26-2elem_slowsys_integrate_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)