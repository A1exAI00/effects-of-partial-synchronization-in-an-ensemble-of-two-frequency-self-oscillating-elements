using Statistics: mean

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

""" Не применять для определения пересечения секущей """
function passing_through_zero_positive_i(arr)
    idx = []
    for i in 1:length(arr)-1
        if (arr[i] ≤ 0) && (arr[i+1] > 0)
            push!(idx, i)
        end
    end
    return idx
end

""" Не применять для определения пересечения секущей """
function passing_through_zero_negative_i(arr)
    idx = []
    for i in 1:length(arr)-1
        if (arr[i] ≥ 0) && (arr[i+1] < 0)
            push!(idx, i)
        end
    end
    return idx
end

""" Улучшеная версия алгоритма определения времени пересечения секущей 
с использованием линейной интерполяции """
function passing_through_zero_positive_t(arr, t_arr)
    t_pass = Float64[]

    for i in 1:length(arr)-1
        if (arr[i] ≤ 0) && (arr[i+1] > 0)
            x₁, x₂ = arr[i], arr[i+1]
            t₁, t₂ = t_arr[i], t_arr[i+1]
            k = (x₂-x₁)/(t₂-t₁)
            b = x₁ - k*t₁
            t₀ = -b/k
            push!(t_pass, t₀)
        end
    end
    return t_pass
end

function calc_avg_period(arr, t_arr)
    t_pass_zero = passing_through_zero_positive_t(arr, t_arr)
    if length(t_pass_zero) < 2 return NaN end
    return mean(diff(t_pass_zero))
end

function calc_avg_freq(arr, t_arr)
    return 1/calc_avg_period(arr, t_arr)
end

""" TODO: заменить индексацию на линейную интеполяцию. 
Хотя в принципе когда u пересекает значение 0.0, значение v остается 
примерно одинаковым, так что можно оставить интексацию """
function calc_avg_amtlitude(arr_u, arr_v)
    idx_pass_zero_v = passing_through_zero_negative_i(arr_v)
    if length(idx_pass_zero_v) < 2 return NaN end
    return mean(arr_u[idx_pass_zero_v])
end

function calc_phase(tₙ::Vector{Float64}, t::Float64)
    if length(tₙ) < 2 
        return NaN
    end

    if t in tₙ
        return 0.0
    end

    t_last::Float64, t_next::Float64 = 0.0, 0.0
    if (t < tₙ[1])
        t_last, t_next = tₙ[1], tₙ[2]
        curr_ω = 1/(t_next-t_last)
        φ = 2π * curr_ω * (t-t_last)
        return rem2pi(φ, RoundNearest)
    end
    if t > tₙ[end]
        t_last, t_next = tₙ[end-1], tₙ[end]
        curr_ω = 1/(t_next-t_last)
        φ = 2π * curr_ω * (t-t_last)
        return rem2pi(φ, RoundNearest)
    end

    for i in 1:(length(tₙ)-1)
        if (t > tₙ[i]) && (t < tₙ[i+1])
            t_last, t_next = tₙ[i], tₙ[i+1]
            break
        end
    end

    curr_ω = 1/(t_next-t_last)
    φ = 2π * curr_ω * (t-t_last)
    return rem2pi(φ, RoundNearest)
end

function calc_phase(arr::Vector{Float64}, t_arr::Vector{Float64}, t::Float64)
    tₙ = passing_through_zero_positive_t(arr, t_arr)
    return calc_phase(tₙ, t)
end

#########################################################################################

function last(arr, size)
    if length(arr) ≤ size
        return arr
    else
        return arr[end-size+1:end]
    end
end