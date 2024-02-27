
using DifferentialEquations

function nonlinearity_1(u, a, b, c)
    return -u*(u+a)*(u-a)*(u+b)*(u-b)*(u+c)*(u-c)
end


function article_model_single_element(du, u, p, t)
    a, b, c, ε, J = p
    x, y = u
    du[1] = nonlinearity_1(x, a, b, c) - y
    du[2] = ε*(x-J)
    return nothing
end

function integrate(u₀, t_span, p)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    # ALG = Tsit5()
    ALG = Rosenbrock23()

    prob = ODEProblem(article_model_single_element, u₀, t_span, p)
    sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)

    return (sol[1, :], sol[2, :])
end