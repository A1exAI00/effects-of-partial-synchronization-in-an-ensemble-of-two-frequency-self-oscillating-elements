#=
Функционал связанный с начальными условиями для численного интегрирования
=#

#########################################################################################

Random.seed!(42)

IC_N = 200 # Количество начальных условий

a_init_x = [range(A₃[1], B₃[1], round(Int, IC_N/2))..., range(D₃[1], C₃[1], round(Int, IC_N/2))...]
a_init_y = f.(a_init_x)
a_init = collect(zip(a_init_x, a_init_y))

c_init_x = [range(A₁[1], B₁[1], round(Int, IC_N/2))..., range(D₁[1], C₁[1], round(Int, IC_N/2))...]
c_init_y = f.(c_init_x)
c_init = collect(zip(c_init_x, c_init_y))

#########################################################################################

function initial_random_phase(is_c)
    if is_c
        return shuffle(c_init)[1]
    end
    return shuffle(a_init)[1]
end

function initial_zero_phase(is_c)
    if is_c
        return (c, 0)
    end
    return (a, 0)
end

function initial_ring_phase(is_c, i)
    if is_c
        return c_init[rem(i, IC_N)]
    end
    return a_init[rem(i, IC_N)]
end

function initial_asinphase(is_c, i)
    i = iseven(i) ? rem(i, IC_N) : rem(i+round(Int, IC_N/2), IC_N)
    if is_c
        return c_init[i]
    end
    return a_init[i]
end

#########################################################################################

function φ_mode_to_init_points(φ_mode, init_pattern)
    if φ_mode == "random"
        init_points = initial_random_phase.(init_pattern)
    elseif φ_mode == "zero"
        init_points = initial_zero_phase.(init_pattern)
    elseif φ_mode == "синфазно"
        init_points = [initial_ring_phase(init_pattern[i], i) for i in eachindex(init_pattern)]
    elseif φ_mode == "противофазно"
        init_points = [initial_asinphase(init_pattern[i], i) for i in eachindex(init_pattern)]
    else
        error("Неверно выбран параметр φ_mode=$φ_mode")
    end
    
end