#=
Цель src-файла: определить функции, необходимые для интегрирования системы
=#

#########################################################################################

function model_single_element(du, u, p, t)
    a, b, c, ε, J = p
    x, y = u
    du[1] = f(x,a,b,c) - y
    du[2] = ε*(x-J)
    return nothing
end

function integrate_single_element(u₀, t_span, J)
    ABSTOL = 1e-7
    RELTOL = 1e-7
    # ALG = Tsit5()
    ALG = Rosenbrock23()

    prob = ODEProblem(model_single_element, u₀, t_span, (a,b,c,ε,J))
    sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)
    return sol
end

#########################################################################################

"""
dU = [duᵢ..., dvᵢ...] \\
U = [uᵢ..., vᵢ...] \\
p = [a, b, c, ε, d, N]
"""
function model_multiple_elements(dU, U, p, t)
    a, b, c, ε, d, N = p
    u = U[1:N]
    v = U[N+1:2N]

    for i in 1:N
        if i == 1
            dU[i] = f(u[i], a, b, c) - v[i] + d*(u[N]-2*u[1]+u[2])
        elseif i == N
            dU[i] = f(u[i], a, b, c) - v[i] + d*(u[N-1]-2*u[N]+u[1])
        else
            dU[i] = f(u[i], a, b, c) - v[i] + d*(u[i-1]-2*u[i]+u[i+1])
        end
        dU[N+i] = ε*(u[i])
    end
    return nothing
end

function integrate_multiple_elements(U₀, t_span, d, N_elements; saveat=nothing)
    ABSTOL = 1e-3
    RELTOL = 1e-3
    # ALG = Tsit5()
    ALG = Rosenbrock23()
    MAXITERS = Int(1e9)

    p = (a, b, c, ε, d, N_elements)

    prob = ODEProblem(model_multiple_elements, U₀, t_span, p)
    if isnothing(saveat)
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS)
    else
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS, saveat=saveat)
    end

    return sol
end

#########################################################################################

# function article1_calc_avg_freq(U₀, p, t_span)
#     sol = integrate_multiple_elements(U₀, p, t_span)
#     sol_t = sol.t
#     uᵢ = [sol[i,:] for i in 1:N_elements]

#     ωᵢ = []
#     for i in 1:N_elements
#         push!(ωᵢ, calc_avg_freq(uᵢ[i], sol_t))
#     end
#     return ωᵢ
# end

# function article1_calc_avg_amplitude(U₀, p, t_span)
#     sol = integrate_multiple_elements(U₀, p, t_span)
#     uᵢ = [sol[i,:] for i in 1:N_elements]
#     vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

#     aᵢ = []
#     for i in 1:N_elements
#         push!(aᵢ, calc_avg_amtlitude(vᵢ[i], uᵢ[i]))
#     end
#     return aᵢ
# end

# function atricle1_calc_final_phase(U₀, p, t_span)
#     sol = integrate_multiple_elements(U₀, p, t_span)
#     sol_t = sol.t
#     uᵢ = [sol[i,:] for i in 1:N_elements]
#     vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

#     φᵢ = []
#     for i in 1:N_elements
#         push!(φᵢ, calc_phase(uᵢ[i], sol_t, sol_t[end]))
#     end
#     return φᵢ
# end