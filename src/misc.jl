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

function passing_through_zero_positive_i(arr)
    idx = []
    for i in 1:length(arr)-1
        if (arr[i] ≤ 0) && (arr[i+1] > 0)
            push!(idx, i)
        end
    end
    return idx
end

function passing_through_zero_negative_i(arr)
    idx = []
    for i in 1:length(arr)-1
        if (arr[i] ≥ 0) && (arr[i+1] < 0)
            push!(idx, i)
        end
    end
    return idx
end

function passing_through_zero_positive_t(arr, t_arr)
    idx_pass_zero = passing_through_zero_positive_i(arr)
    t_pass_zero = t_arr[idx_pass_zero]
    return t_pass_zero
end

function calc_period(arr, t_arr)
    idx_pass_zero = passing_through_zero_positive_i(arr)
    if length(idx_pass_zero) > 10
        idx_middle = round(Int, length(idx_pass_zero)/2)
    else
        error("not enough time to calc period")
    end
    t_pass_zero = t_arr[idx_pass_zero]
    t_pass_zero = t_pass_zero[idx_pass_zero .> idx_middle]
    return mean(diff(t_pass_zero))
end

function calc_amtlitude(arr_u, arr_v)
    idx_pass_zero_v = passing_through_zero_negative_i(arr_v)
    if length(idx_pass_zero_v) > 10
        idx_middle = round(Int, length(idx_pass_zero_v)/2)
    else
        error("not enough time to calc period")
    end
    idx_pass_zero_v = idx_pass_zero_v[idx_pass_zero_v .> idx_middle]
    return mean(arr_u[idx_pass_zero_v])
end

#########################################################################################

# test calc_period(arr, t_arr)
# x = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]
# t = eachindex(x)
# println(calc_period(x,t))