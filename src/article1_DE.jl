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

function integrate_multiple_elements(U₀, t_span, d, N_elements; saveat=nothing, alg=nothing)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    MAXITERS = Int(1e9)
    # ALG = Tsit5()
    # ALG = Rosenbrock23()
    # ALG = lsoda() 
    if isnothing(alg)
        ALG = FBDF()
    else
        ALG = alg
    end

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

function moded_model_multiple_elements(dU, U, p, t)
    a, b, c, ε, d, N, J = p
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
        dU[N+i] = ε*(u[i]-J)
    end
    return nothing
end

function moded_integrate_multiple_elements(U₀, t_span, d, N_elements, J; saveat=nothing, alg=nothing)
    ABSTOL = 1e-3
    RELTOL = 1e-3
    MAXITERS = Int(1e9)
    # ALG = Tsit5()
    # ALG = Rosenbrock23()
    # ALG = lsoda()
    if isnothing(alg)
        ALG = FBDF()
    else
        ALG = alg
    end

    p = (a, b, c, ε, d, N_elements, J)

    prob = ODEProblem(moded_model_multiple_elements, U₀, t_span, p)
    if isnothing(saveat)
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS)
    else
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS, saveat=saveat)
    end

    return sol
end

#########################################################################################

function moded_model_two_elements(dU, U, p, t)
    a, b, c, ε, d, J = p
    u₁, u₂, v₁, v₂ = U

    dU[1] = f(u₁, a, b, c) - v₁ + d*(u₂-u₁)
    dU[2] = f(u₂, a, b, c) - v₂ + d*(u₁-u₂)
    dU[3] = ε*(u₁-J)
    dU[4] = ε*(u₂-J)

    return nothing
end

function moded_integrate_two_elements(U₀, t_span, d, J; saveat=nothing, alg=nothing)
    ABSTOL = 1e-3
    RELTOL = 1e-3
    MAXITERS = Int(1e9)
    # ALG = Tsit5()
    # ALG = Rosenbrock23()
    # ALG = lsoda()
    if isnothing(alg)
        ALG = FBDF()
    else
        ALG = alg
    end

    p = (a, b, c, ε, d, J)

    prob = ODEProblem(moded_model_two_elements, U₀, t_span, p)
    if isnothing(saveat)
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS)
    else
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS, saveat=saveat)
    end

    return sol
end