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

# x range for nonlinearity
x_from_c_offset = 0.05
x_start = -article1.c - x_from_c_offset
x_end = article1.c - x_from_c_offset
Δx = 1e-5
x_N = round(Int,(x_end - x_start)/Δx)
x_range = range(x_start, x_end, x_N)

# nonlinearity
f_res = article1.f.(x_range)
x_extr = [article1.D₃[1], article1.C₂[1], article1.D₁[1]]
y_extr = [article1.D₃[2], article1.C₂[2], article1.D₁[2]]

t_start, t_end = 0.0, 2e4
t_span = (t_start, t_end)

#########################################################################################

J₁ = 0
J₂ = 0.165
J₃ = 0.173
J₄ = 0.176
J₅ = 0.4
J₆ = 0.6

#########################################################################################

sol_a_J₁ = article1.integrate_single_element([ article1.a,0], t_span, J₁)
sol_b_J₁ = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J₁)
sol_c_J₁ = article1.integrate_single_element([ article1.c,0], t_span, J₁)
sol_a_x_J₁, sol_a_y_J₁ = sol_a_J₁[1, :], sol_a_J₁[2, :]
sol_b_x_J₁, sol_b_y_J₁ = sol_b_J₁[1, :], sol_b_J₁[2, :]
sol_c_x_J₁, sol_c_y_J₁ = sol_c_J₁[1, :], sol_c_J₁[2, :]

sol_a_J₂ = article1.integrate_single_element([ article1.a,0], t_span, J₂)
sol_b_J₂ = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J₂)
sol_c_J₂ = article1.integrate_single_element([ article1.c,0], t_span, J₂)
sol_a_x_J₂, sol_a_y_J₂ = sol_a_J₂[1, :], sol_a_J₂[2, :]
sol_b_x_J₂, sol_b_y_J₂ = sol_b_J₂[1, :], sol_b_J₂[2, :]
sol_c_x_J₂, sol_c_y_J₂ = sol_c_J₂[1, :], sol_c_J₂[2, :]

sol_a_J₃ = article1.integrate_single_element([ article1.a,0], t_span, J₃)
sol_b_J₃ = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J₃)
sol_c_J₃ = article1.integrate_single_element([ article1.c,0], t_span, J₃)
sol_a_x_J₃, sol_a_y_J₃ = sol_a_J₃[1, :], sol_a_J₃[2, :]
sol_b_x_J₃, sol_b_y_J₃ = sol_b_J₃[1, :], sol_b_J₃[2, :]
sol_c_x_J₃, sol_c_y_J₃ = sol_c_J₃[1, :], sol_c_J₃[2, :]

sol_a_J₄ = article1.integrate_single_element([ article1.a,0], t_span, J₄)
sol_b_J₄ = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J₄)
sol_c_J₄ = article1.integrate_single_element([ article1.c,0], t_span, J₄)
sol_a_x_J₄, sol_a_y_J₄ = sol_a_J₄[1, :], sol_a_J₄[2, :]
sol_b_x_J₄, sol_b_y_J₄ = sol_b_J₄[1, :], sol_b_J₄[2, :]
sol_c_x_J₄, sol_c_y_J₄ = sol_c_J₄[1, :], sol_c_J₄[2, :]

sol_a_J₅ = article1.integrate_single_element([ article1.a,0], t_span, J₅)
sol_b_J₅ = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J₅)
sol_c_J₅ = article1.integrate_single_element([ article1.c,0], t_span, J₅)
sol_a_x_J₅, sol_a_y_J₅ = sol_a_J₅[1, :], sol_a_J₅[2, :]
sol_b_x_J₅, sol_b_y_J₅ = sol_b_J₅[1, :], sol_b_J₅[2, :]
sol_c_x_J₅, sol_c_y_J₅ = sol_c_J₅[1, :], sol_c_J₅[2, :]

sol_a_J₆ = article1.integrate_single_element([ article1.a,0], t_span, J₆)
sol_b_J₆ = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J₆)
sol_c_J₆ = article1.integrate_single_element([ article1.c,0], t_span, J₆)
sol_a_x_J₆, sol_a_y_J₆ = sol_a_J₆[1, :], sol_a_J₆[2, :]
sol_b_x_J₆, sol_b_y_J₆ = sol_b_J₆[1, :], sol_b_J₆[2, :]
sol_c_x_J₆, sol_c_y_J₆ = sol_c_J₆[1, :], sol_c_J₆[2, :]

#########################################################################################

fig = Figure(size=(700, 1050))
ax1 = beautiful_Axis(fig[1, 1]; title="J=$(J₁)", xlabel="u", ylabel="v")
ax2 = beautiful_Axis(fig[1, 2]; title="J=$(J₂)", xlabel="u", ylabel="v")
ax3 = beautiful_Axis(fig[2, 1]; title="J=$(J₃)", xlabel="u", ylabel="v")
ax4 = beautiful_Axis(fig[2, 2]; title="J=$(J₄)", xlabel="u", ylabel="v")
ax5 = beautiful_Axis(fig[3, 1]; title="J=$(J₅)", xlabel="u", ylabel="v")
ax6 = beautiful_Axis(fig[3, 2]; title="J=$(J₆)", xlabel="u", ylabel="v")

vlines!.((ax1, ax2, ax3, ax4, ax5, ax6), 0.0, color=:black)
hlines!.((ax1, ax2, ax3, ax4, ax5, ax6), 0.0, color=:black)

vlines!(ax1, J₁, color=:blue, linestyle=:dash)
vlines!(ax2, J₂, color=:blue, linestyle=:dash)
vlines!(ax3, J₃, color=:blue, linestyle=:dash)
vlines!(ax4, J₄, color=:blue, linestyle=:dash)
vlines!(ax5, J₅, color=:blue, linestyle=:dash)
vlines!(ax6, J₆, color=:blue, linestyle=:dash)

lines!(ax1, x_range, f_res, label="f(u)", color=:green)
lines!(ax2, x_range, f_res, label="f(u)", color=:green)
lines!(ax3, x_range, f_res, label="f(u)", color=:green)
lines!(ax4, x_range, f_res, label="f(u)", color=:green)
lines!(ax5, x_range, f_res, label="f(u)", color=:green)
lines!(ax6, x_range, f_res, label="f(u)", color=:green)


lines!(ax1, sol_a_x_J₁[round(Int, end/2):end], sol_a_y_J₁[round(Int, end/2):end], color=:blue)
lines!(ax1, sol_b_x_J₁[round(Int, end/2):end], sol_b_y_J₁[round(Int, end/2):end], color=:red)
lines!(ax1, sol_c_x_J₁[round(Int, end/2):end], sol_c_y_J₁[round(Int, end/2):end], color=:blue)

lines!(ax2, sol_a_x_J₂[round(Int, end/2):end], sol_a_y_J₂[round(Int, end/2):end], color=:blue)
lines!(ax2, sol_b_x_J₂[round(Int, end/2):end], sol_b_y_J₂[round(Int, end/2):end], color=:red)
lines!(ax2, sol_c_x_J₂[round(Int, end/2):end], sol_c_y_J₂[round(Int, end/2):end], color=:blue)

lines!(ax3, sol_a_x_J₃[round(Int, end/2):end], sol_a_y_J₃[round(Int, end/2):end], color=:blue)
lines!(ax3, sol_b_x_J₃[round(Int, end/2):end], sol_b_y_J₃[round(Int, end/2):end], color=:red)
lines!(ax3, sol_c_x_J₃[round(Int, end/2):end], sol_c_y_J₃[round(Int, end/2):end], color=:blue)

lines!(ax4, sol_a_x_J₄[round(Int, end/2):end], sol_a_y_J₄[round(Int, end/2):end], color=:blue)
lines!(ax4, sol_b_x_J₄[round(Int, end/2):end], sol_b_y_J₄[round(Int, end/2):end], color=:red)
lines!(ax4, sol_c_x_J₄[round(Int, end/2):end], sol_c_y_J₄[round(Int, end/2):end], color=:blue)

lines!(ax5, sol_a_x_J₅[round(Int, end/2):end], sol_a_y_J₅[round(Int, end/2):end], color=:blue)
lines!(ax5, sol_b_x_J₅[round(Int, end/2):end], sol_b_y_J₅[round(Int, end/2):end], color=:red)
lines!(ax5, sol_c_x_J₅[round(Int, end/2):end], sol_c_y_J₅[round(Int, end/2):end], color=:blue)

lines!(ax6, sol_a_x_J₆[round(Int, end/2):end], sol_a_y_J₆[round(Int, end/2):end], color=:blue)
lines!(ax6, sol_b_x_J₆[round(Int, end/2):end], sol_b_y_J₆[round(Int, end/2):end], color=:red)
lines!(ax6, sol_c_x_J₆[round(Int, end/2):end], sol_c_y_J₆[round(Int, end/2):end], color=:blue)

# axislegend.((ax1, ax2, ax3, ax4, ax5, ax6), position=:rt) # (l, r, c), (b, t, c)

savepath = plotsdir("tmp", "08-Nelem_PS_LaTeX_$(time_ns()).png")
save(savepath, fig, px_per_unit=1.5)
