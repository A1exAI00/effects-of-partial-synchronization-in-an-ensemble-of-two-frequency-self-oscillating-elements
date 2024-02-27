using DrWatson
@quickactivate "semester8"

using CairoMakie

include(srcdir("article_phase_space.jl"))

#########################################################################################

function find_extremum_index(arr)
    extr_i = []
    for i in 2:length(arr)-1
        if ((arr[i] < arr[i-1]) && (arr[i] < arr[i+1])) || ((arr[i] > arr[i-1]) && (arr[i] > arr[i+1]))
            push!(extr_i, i)
        end
    end
    return extr_i
end

#########################################################################################

function plot_phase_space(J; plot_in=true, plot_mid=true, plot_out=true)
    p = (a, b, c, ε, J)
    U₀_lower_1 = [u₀_lower_1, v₀_lower_1]
    U₀_lower_2 = [u₀_lower_2, v₀_lower_2]
    U₀_upper_1 = [u₀_upper_1, v₀_upper_1]
    U₀_upper_2 = [u₀_upper_2, v₀_upper_2]
    t_span = (t_start, t_end)

    

    println("start solve, $(J)")
    x_sol_lower_1, y_sol_lower_1 = integrate(U₀_lower_1, t_span, p)
    x_sol_lower_2, y_sol_lower_2 = integrate(U₀_lower_2, t_span, p)
    x_sol_upper_1, y_sol_upper_1 = integrate(U₀_upper_1, t_span, p)
    x_sol_upper_2, y_sol_upper_2 = integrate(U₀_upper_2, t_span, p)
    x_sol_back_1, y_sol_back_1 = integrate([b,0], (-1).*t_span, p)
    x_sol_back_2, y_sol_back_2 = integrate([-b,0], (-1).*t_span, p)
    # println("end solve")

    

    fig = Figure(size=(1000, 1000))
    ax = Axis(fig[1, 1], 
        title="Фазовое пространство одного элемента: J=$(J)", 
        xlabel="uⱼ", 
        ylabel="vⱼ",
        xminorticksvisible = true, 
        xminorgridvisible = true, 
        yminorticksvisible = true, 
        yminorgridvisible = true, 
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10)
    )

    vlines!(ax, 0.0, color=:black)
    hlines!(ax, 0.0, color=:black)

    lines!(ax, x_range, y_nonlinearity, label="Нелинейность", color=:green)
    vlines!(ax, J, color=:blue, linestyle=:dash)
    scatter!(ax, x_extr, y_extr, label="Экстремумы", color=:green)

    if plot_in
        lines!(ax, x_sol_lower_1, y_sol_lower_1, color=:blue)
        lines!(ax, x_sol_lower_2, y_sol_lower_2, color=:blue)
    end
    if plot_mid
        lines!(ax, x_sol_back_1, y_sol_back_1, color=:red)
        lines!(ax, x_sol_back_2, y_sol_back_2, color=:red)
    end
    if plot_out
        lines!(ax, x_sol_upper_1, y_sol_upper_1, color=:blue)
        lines!(ax, x_sol_upper_2, y_sol_upper_2, color=:blue)
    end

    axislegend(ax, position=:rt) # (l, r, c), (b, t, c)
    save(plotsdir("01-article_phase_space_$(time_ns()).png"), fig, px_per_unit=2)
end

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001

# x range for nonlinearity
x_from_c_offset = 0.05
x_start, x_end = -c-x_from_c_offset, c+x_from_c_offset
Δx = 1e-5
x_N = round(Int,(x_end - x_start)/Δx)
x_range = range(x_start, x_end, x_N)

# calc nonlinearity
nonlinearity_(x) = nonlinearity_1(x, a, b, c)
y_nonlinearity = nonlinearity_.(x_range)

# all extrema of nonlinearity
extr_index = find_extremum_index(y_nonlinearity)
x_extr = x_range[extr_index]
y_extr = y_nonlinearity[extr_index]

# extrema for positive x
first_extr = (x_extr[4], y_extr[4])
second_extr = (x_extr[5], y_extr[5])
third_extr = (x_extr[6], y_extr[6])

#########################################################################################

Js = [
    0, 
    first_extr[1] - 0.01,
    first_extr[1] - 0.001,
    first_extr[1] - 0.0001,
    first_extr[1] - 0.00001,
    first_extr[1],
    first_extr[1] + 0.0001,
    first_extr[1] + 0.001,
    first_extr[1] + 0.01,

    second_extr[1] - 0.01,
    second_extr[1] - 0.001,
    second_extr[1] - 0.0001,
    second_extr[1] - 0.00001,
    second_extr[1],
    second_extr[1] + 0.0001,
    second_extr[1] + 0.001,
    second_extr[1] + 0.01,

    third_extr[1] - 0.01,
    third_extr[1] - 0.001,
    third_extr[1] - 0.0001,
    third_extr[1] - 0.00001,
    third_extr[1],
    third_extr[1] + 0.0001,
    third_extr[1] + 0.001,
    third_extr[1] + 0.01,
]

u₀_lower_1, v₀_lower_1 = a, 0
u₀_lower_2, v₀_lower_2 = -a, 0
u₀_upper_1, v₀_upper_1 = c, 0
u₀_upper_2, v₀_upper_2 = -c, 0

t_start, t_end = 0, 10000

#########################################################################################

for (i,J) in enumerate(Js)
    if i ≤ 20
        plot_phase_space(J)
    else
        plot_phase_space(J; plot_mid=false)
    end
end