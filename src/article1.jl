
using DifferentialEquations
include("misc.jl")

#########################################################################################

function nonlinearity_1(u, a, b, c)
    return -u*(u+a)*(u-a)*(u+b)*(u-b)*(u+c)*(u-c)
end

function ∂F₁_∂U_at_eq(J)
    a², b², c² = a^2, b^2, c^2
    return a²*b²*c² - 3*(a²*b²+a²*c²+b²*c²)*J^2 + 5*(a²+b²+c²)*J^4 - 7*J^6
end

function ∂²F₁_∂U²_at_eq(J)
    a², b², c² = a^2, b^2, c^2
    return -6*(a²*b²+a²*c²+b²*c²)*J + 20*(a²+b²+c²)*J^3 - 42*J^5
end

function ∂³F₁_∂U³_at_eq(J)
    a², b², c² = a^2, b^2, c^2
    return -6*(a²*b²+a²*c²+b²*c²) + 60*(a²+b²+c²)*J^2 - 210*J^4
end

Reλ₁(γ, ε) = (γ^2-4ε>0) ? 0.5*(γ+sqrt(γ^2-4ε)) : 0.5*γ
Reλ₂(γ, ε) = (γ^2-4ε>0) ? 0.5*(γ-sqrt(γ^2-4ε)) : 0.5*γ
Imλ₁(γ, ε) = (γ^2-4ε>0) ? 0 : 0.5*sqrt(4ε-γ^2)
Imλ₂(γ, ε) = (γ^2-4ε>0) ? 0 : -0.5*sqrt(4ε-γ^2)

function first_Lyapunov_quantity(J)
    
    A = ∂F₁_∂U_at_eq(J)
    B = -1
    C = ε
    a₂₀ = ∂²F₁_∂U²_at_eq(J)
    a₃₀ = ∂³F₁_∂U³_at_eq(J)

    # p = - A
    q = -B*C

    return -π/(4*B*q*sqrt(q)) * (-2*A*B*(a₂₀^2) - (A^2 + B*C)*3*(-B*a₃₀))
end

#########################################################################################

function article1_model_single_element(du, u, p, t)
    a, b, c, ε, J = p
    x, y = u
    du[1] = nonlinearity_1(x, a, b, c) - y
    du[2] = ε*(x-J)
    return nothing
end

function integrate(u₀, t_span, p)
    ABSTOL = 1e-7
    RELTOL = 1e-7
    # ALG = Tsit5()
    ALG = Rosenbrock23()

    prob = ODEProblem(article1_model_single_element, u₀, t_span, p)
    sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)

    return (sol[1, :], sol[2, :])
end

#########################################################################################

"""
dU = [duᵢ..., dvᵢ...]
U = [uᵢ..., vᵢ...]
p = [a, b, c, ε, d, N]
"""
function article1_model_multiple_elements(dU, U, p, t)
    a, b, c, ε, d, N = p
    u = U[1:N]
    v = U[N+1:2N]

    for i in 1:N
        if i == 1
            dU[i] = nonlinearity_1(u[i], a, b, c) - v[i] + d*(u[i]-2*u[i]+u[i+1])
        elseif i == N
            dU[i] = nonlinearity_1(u[i], a, b, c) - v[i] + d*(u[i-1]-2*u[i]+u[i])
        else
            dU[i] = nonlinearity_1(u[i], a, b, c) - v[i] + d*(u[i-1]-2*u[i]+u[i+1])
        end
        dU[N+i] = ε*(u[i])
    end
    return nothing
end

function integrate_chain(U₀, p, t_span; saveat=nothing)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    # ALG = Tsit5()
    ALG = Rosenbrock23()

    prob = ODEProblem(article1_model_multiple_elements, U₀, t_span, p)
    if isnothing(saveat)
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)
    else
        sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, saveat=saveat)
    end

    return sol
end

function article1_calc_avg_freq(U₀, p, t_span)
    sol = integrate_chain(U₀, p, t_span)
    sol_t = sol.t
    uᵢ = [sol[i,:] for i in 1:N_elements]

    ωᵢ = []
    for i in 1:N_elements
        push!(ωᵢ, calc_avg_freq(uᵢ[i], sol_t))
    end
    return ωᵢ
end

function article1_calc_avg_amplitude(U₀, p, t_span)
    sol = integrate_chain(U₀, p, t_span)
    uᵢ = [sol[i,:] for i in 1:N_elements]
    vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

    aᵢ = []
    for i in 1:N_elements
        push!(aᵢ, calc_avg_amtlitude(vᵢ[i], uᵢ[i]))
    end
    return aᵢ
end

function atricle1_calc_final_phase(U₀, p, t_span)
    sol = integrate_chain(U₀, p, t_span)
    sol_t = sol.t
    uᵢ = [sol[i,:] for i in 1:N_elements]
    vᵢ = [sol[N_elements+i,:] for i in 1:N_elements]

    φᵢ = []
    for i in 1:N_elements
        push!(φᵢ, calc_phase(uᵢ[i], sol_t, sol_t[end]))
    end
    return φᵢ
end