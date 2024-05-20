#=

=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie
using NLsolve

include(srcdir("article1_module.jl"))
include(srcdir("misc.jl"))
include(srcdir("plotting_tools.jl"))

using .article1

#########################################################################################

d_start, d_end, d_N = 0.0, 0.4, 500
J_start, J_end, J_N = 0.0, article1.c, 100

#########################################################################################

f(x) = article1.f(x)
dfdx(x) = article1.f_deriv(x)

function opt_Imλ12!(F,x)
    F[1] = abs(dfdx(x[1])) - 2*sqrt(article1.ε)
end

function gen_opt_Reλ34(d)
    function opt_Reλ34!(F,x)
        F[1] = dfdx(x[1]) - 2*d
    end
    return opt_Reλ34!
end

function gen_opt_Imλ34(d)
    function opt_Imλ34!(F,x)
        F[1] = abs(dfdx(x[1])-2*d) - 2*sqrt(article1.ε)
    end
    return opt_Imλ34!
end

#########################################################################################

d_range = range(d_start, d_end, d_N)
J_range = range(J_start, J_end, J_N)

J_extr = [article1.D₃[1], article1.C₂[1], article1.D₁[1]]

#########################################################################################

J_Imλ12 = []
J_Reλ34 = []
J_Imλ34 = []

for (i,d) in enumerate(d_range)
    J_Imλ12_curr = zeros(6)
    J_Imλ12_inits = [J_extr[1] - 1e-1, J_extr[1] + 1e-1,
                    J_extr[2] - 1e-1, J_extr[2] + 1e-1,
                    J_extr[3] - 1e-2, J_extr[3] + 1e-2,
    ]
    J_Imλ12_curr[1] = nlsolve(opt_Imλ12!, [J_Imλ12_inits[1]], autodiff=:forward).zero[1]
    J_Imλ12_curr[2] = nlsolve(opt_Imλ12!, [J_Imλ12_inits[2]], autodiff=:forward).zero[1]
    J_Imλ12_curr[3] = nlsolve(opt_Imλ12!, [J_Imλ12_inits[3]], autodiff=:forward).zero[1]
    J_Imλ12_curr[4] = nlsolve(opt_Imλ12!, [J_Imλ12_inits[4]], autodiff=:forward).zero[1]
    J_Imλ12_curr[5] = nlsolve(opt_Imλ12!, [J_Imλ12_inits[5]], autodiff=:forward).zero[1]
    J_Imλ12_curr[6] = nlsolve(opt_Imλ12!, [J_Imλ12_inits[6]], autodiff=:forward).zero[1]
    push!(J_Imλ12, J_Imλ12_curr)

    J_Reλ34_prev = i==1 ? J_extr : J_Reλ34[end]
    J_Reλ34_curr = zeros(3)
    J_Reλ34_curr[1] = nlsolve(gen_opt_Reλ34(d), [J_Reλ34_prev[1]], autodiff=:forward).zero[1]
    J_Reλ34_curr[2] = nlsolve(gen_opt_Reλ34(d), [J_Reλ34_prev[2]], autodiff=:forward).zero[1]
    J_Reλ34_curr[3] = nlsolve(gen_opt_Reλ34(d), [J_Reλ34_prev[3]], autodiff=:forward).zero[1]
    push!(J_Reλ34, J_Reλ34_curr)

    J_Imλ34_curr = zeros(6)
    J_Imλ34_inits = (i==1) ? J_Imλ12_curr : J_Imλ34[end]
    J_Imλ34_curr[1] = nlsolve(gen_opt_Imλ34(d), [J_Imλ34_inits[1]], autodiff=:forward).zero[1]
    J_Imλ34_curr[2] = nlsolve(gen_opt_Imλ34(d), [J_Imλ34_inits[2]], autodiff=:forward).zero[1]
    J_Imλ34_curr[3] = nlsolve(gen_opt_Imλ34(d), [J_Imλ34_inits[3]], autodiff=:forward).zero[1]
    J_Imλ34_curr[4] = nlsolve(gen_opt_Imλ34(d), [J_Imλ34_inits[4]], autodiff=:forward).zero[1]
    J_Imλ34_curr[5] = nlsolve(gen_opt_Imλ34(d), [J_Imλ34_inits[5]], autodiff=:forward).zero[1]
    J_Imλ34_curr[6] = nlsolve(gen_opt_Imλ34(d), [J_Imλ34_inits[6]], autodiff=:forward).zero[1]
    push!(J_Imλ34, J_Imλ34_curr)
end

#########################################################################################

fig = Figure(size=(700, 500))
ax = Axis(fig[1, 1]; title="", xlabel="J", ylabel="d")

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

lines!(ax, J_extr[1].+0*d_range, d_range, color=:blue, label="Reλ₁₂=0")
lines!(ax, J_extr[2].+0*d_range, d_range, color=:blue)
lines!(ax, J_extr[3].+0*d_range, d_range, color=:blue)

# lines!(ax, d_range, [J_Imλ12[j][1] for j in eachindex(d_range)], color=:blue, linestyle=:dash, label="Imλ₁₂=0")
# lines!(ax, d_range, [J_Imλ12[j][2] for j in eachindex(d_range)], color=:blue, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ12[j][3] for j in eachindex(d_range)], color=:blue, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ12[j][4] for j in eachindex(d_range)], color=:blue, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ12[j][5] for j in eachindex(d_range)], color=:blue, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ12[j][6] for j in eachindex(d_range)], color=:blue, linestyle=:dash)

lines!(ax, [J_Reλ34[j][1] for j in eachindex(d_range)], d_range, color=:red, label="Reλ₃₄=0")
lines!(ax, [J_Reλ34[j][2] for j in eachindex(d_range)], d_range, color=:red)
lines!(ax, [J_Reλ34[j][3] for j in eachindex(d_range)], d_range, color=:red)

# lines!(ax, d_range, [J_Imλ34[j][1] for j in eachindex(d_range)], color=:red, linestyle=:dash, label="Imλ₃₄=0")
# lines!(ax, d_range, [J_Imλ34[j][2] for j in eachindex(d_range)], color=:red, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ34[j][3] for j in eachindex(d_range)], color=:red, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ34[j][4] for j in eachindex(d_range)], color=:red, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ34[j][5] for j in eachindex(d_range)], color=:red, linestyle=:dash)
# lines!(ax, d_range, [J_Imλ34[j][6] for j in eachindex(d_range)], color=:red, linestyle=:dash)

limits!(ax, -0.05, nothing, -0.01, nothing)
axislegend(ax, position=:rb) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "27-2elem_(J,d)_plane_lambdas_$(time_ns()).png")
save(savepath, fig, px_per_unit=2)