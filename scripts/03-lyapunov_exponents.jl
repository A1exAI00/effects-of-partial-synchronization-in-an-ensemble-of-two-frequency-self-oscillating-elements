#=
Цель программы: изобразить график зависимости ляпуновских показателей 
в зависимости от значения параметра J.

Для определения устойчивости СР при значениях J, при которых наблюдается 
бифуркация Андронова-Хопфа, линеаризованной системы недостаточно, и систему надо 
раскладывать в окрестности СР в ряд Тейлора с учётом членов высших порядков (до ненулевого).

Аналитически нормальная форма бифуркации Андронова-Хопфа выглядит следующим образом:

du/dt = α⋅u - β⋅v + L⋅(u²+v²)⋅u
dv/dt = β⋅u + α⋅v + L⋅(u²+v²)⋅v

где L - первая ляпуновская величина.

Для приведения системы к нормальной форме необходимо:
1. Разложить систему в ряд Тейлора до 3-гочлена; 
2. Применить линейное преобразование, которое приведет линейную часть к жордановой форме;
3. Применить нелинейное преобразование, после которого пропадут квадратичные члены

Получить выражение для первой ляпуновской челичины у меня не получилось, из-за того, 
что я остановился на приведении линеаризованной системы к жордановой форме.

Также, насколько я понял, ляпуновскую величину можно получить из реализации.
Изначально я считал, что lyapunov exponents это то же самое, что значение 
первой ляпуновской величины, но кажется это оказалось не так.
=#

#########################################################################################

using DrWatson
@quickactivate "semester8"

using CairoMakie
using DynamicalSystems

include(srcdir("article1.jl"))
include(srcdir("misc.jl"))

#########################################################################################

function calc_lyapunov(a, b, c, ε, J)
    p = (a, b, c, ε, J)
    U₀ = [J, nonlinearity_1(J, a, b, c)]
    ds = CoupledODEs(article_model_single_element, U₀, p)
    return lyapunovspectrum(ds, 10000)
end

#########################################################################################

a, b, c, ε = 0.32, 0.79, 1.166, 0.001

J_start, J_end, J_N = 0.0, c+0.05, 1000
J_range = range(J_start, J_end, J_N)

#########################################################################################

f_ = nonlinearity_1.(J_range, a, b, c)
extr_index = find_extremum_index(f_)
J_extr = J_range[extr_index]
f_extr = f_[extr_index]

# println(calc_lyapunov(a, b, c, ε, 1.0))
lyap_exp = calc_lyapunov.(a, b, c, ε, J_range)

lyap_exp_1 = [lyap_exp[i][1] for i in eachindex(lyap_exp)]
lyap_exp_2 = [lyap_exp[i][2] for i in eachindex(lyap_exp)]

#########################################################################################

fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1], 
    title="Ляпуновские экспоненты от J", 
    xlabel="J", 
    ylabel="Lyapunov exponents",
	xminorticksvisible = true, 
	xminorgridvisible = true, 
	yminorticksvisible = true, 
	yminorgridvisible = true, 
	xminorticks = IntervalsBetween(10),
	yminorticks = IntervalsBetween(10)
)

vlines!(ax, 0.0, color=:black)
hlines!(ax, 0.0, color=:black)

vlines!.(ax, J_extr, color=:green, linestyle=:dash)

lines!(ax, J_range, lyap_exp_1, label="1")
lines!(ax, J_range, lyap_exp_2, label="2")

axislegend(ax, position=:rb) # (l, r, c), (b, t, c)
save(plotsdir("03-lyapunov_exponents_$(time_ns()).png"), fig, px_per_unit=2)