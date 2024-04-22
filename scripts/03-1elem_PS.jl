#=
Цель программы: изобразить фазовое пространство системы.

Для кадого значения J отрисовывается фазовый портрет. На нем изображены 
4 фазовые траектории, интегрированные в прямом времени, выходящие из точек 
(a,0), (-a,0), (c,0), (-c,0), и две фазовые траектории, интегрированные в 
обратном времени, выходящие из точек (b,0), (-b,0).

Эти фазовые траектории наматываются на устойчивые и неустойчивые предельные циклы.

TODO: начиналь отрисовку фазовых траекторий не сразу, а когда придут к устойчивым стуктурам
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

savedir = joinpath("tmp", "03-1elem_PS")
if !(ispath(savedir)) mkpath(plotsdir(savedir)) end

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

Js = [
    0, 
    x_extr[1] - 0.01,
    x_extr[1] - 0.001,
    x_extr[1] - 0.0001,
    x_extr[1] - 0.00001,
    x_extr[1],
    x_extr[1] + 0.0001,
    x_extr[1] + 0.001,
    x_extr[1] + 0.01,

    x_extr[2] - 0.01,
    x_extr[2] - 0.001,
    x_extr[2] - 0.0001,
    x_extr[2] - 0.00001,
    x_extr[2],
    x_extr[2] + 0.0001,
    x_extr[2] + 0.001,
    x_extr[2] + 0.01,

    x_extr[3] - 0.01,
    x_extr[3] - 0.001,
    x_extr[3] - 0.0001,
    x_extr[3] - 0.00001,
    x_extr[3],
    x_extr[3] + 0.0001,
    x_extr[3] + 0.001,
    x_extr[3] + 0.01,
]

#########################################################################################

for (i,J) in enumerate(Js)
    plot_a_limit_cycle = true
    plot_b_limit_cycle = true
    plot_c_limit_cycle = true

    println("start solve, $(J)")
    sol_a_1 = article1.integrate_single_element([ article1.a,0], t_span, J)
    sol_a_2 = article1.integrate_single_element([-article1.a,0], t_span, J)
    sol_b_1 = article1.integrate_single_element([ article1.b,0], (-1).*t_span, J)
    sol_b_2 = article1.integrate_single_element([-article1.b,0], (-1).*t_span, J)
    sol_c_1 = article1.integrate_single_element([ article1.c,0], t_span, J)
    sol_c_2 = article1.integrate_single_element([-article1.c,0], t_span, J)

    sol_a_1_x, sol_a_1_y = sol_a_1[1, :], sol_a_1[2, :]
    sol_a_2_x, sol_a_2_y = sol_a_2[1, :], sol_a_2[2, :]
    sol_b_1_x, sol_b_1_y = sol_b_1[1, :], sol_b_1[2, :]
    sol_b_2_x, sol_b_2_y = sol_b_2[1, :], sol_b_2[2, :]
    sol_c_1_x, sol_c_1_y = sol_c_1[1, :], sol_c_1[2, :]
    sol_c_2_x, sol_c_2_y = sol_c_2[1, :], sol_c_2[2, :]

    fig = Figure(size=(1000, 1000))
    ax = beautiful_Axis(fig[1, 1]; 
        title="Фазовое пространство одного элемента: J=$(J)", 
        xlabel="uⱼ", ylabel="vⱼ",
    )

    vlines!(ax, 0.0, color=:black)
    hlines!(ax, 0.0, color=:black)

    lines!(ax, x_range, f_res, label="Нелинейность", color=:green)
    vlines!(ax, J, color=:blue, linestyle=:dash)
    scatter!(ax, x_extr, y_extr, label="Экстремумы", color=:green)

    if plot_a_limit_cycle
        lines!(ax, sol_a_1_x[round(Int, end/2):end], sol_a_1_y[round(Int, end/2):end], color=:blue)
        lines!(ax, sol_a_2_x[round(Int, end/2):end], sol_a_2_y[round(Int, end/2):end], color=:blue)
    end
    if plot_b_limit_cycle
        lines!(ax, sol_b_1_x[round(Int, end/2):end], sol_b_1_y[round(Int, end/2):end], color=:red)
        lines!(ax, sol_b_2_x[round(Int, end/2):end], sol_b_2_y[round(Int, end/2):end], color=:red)
    end
    if plot_c_limit_cycle
        lines!(ax, sol_c_1_x[round(Int, end/2):end], sol_c_1_y[round(Int, end/2):end], color=:blue)
        lines!(ax, sol_c_2_x[round(Int, end/2):end], sol_c_2_y[round(Int, end/2):end], color=:blue)
    end

    axislegend(ax, position=:rt) # (l, r, c), (b, t, c)

    savepath = plotsdir(savedir, "03-1elem_PS_$(time_ns()).png")
    save(savepath, fig, px_per_unit=2)
end