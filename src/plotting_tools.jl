#=
Цель src-файла: определить функции для графиков
=#

function beautiful_Axis(figpart; title=nothing, xlabel=nothing, ylabel=nothing)
    return Axis(figpart, 
        title=title, 
        xlabel=xlabel, 
        ylabel=ylabel,
        xminorticksvisible = true, 
        xminorgridvisible = true, 
        yminorticksvisible = true, 
        yminorgridvisible = true, 
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10)
    )
end