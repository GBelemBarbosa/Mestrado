include("variables_logistic.jl")
using BenchmarkTools

function ∇f(θ:: Vector{Float64}; X=X, y=y, λ=λ)
    aux=X'y
    e=exp.(-θ.*aux)
    return θ.*(λ.+aux.^2 .*e./(1 .+e))
end

function ∇g(θ:: Vector{Float64}; X=X, y=y, λ=λ)
    aux=X'y
    return θ.*(λ.+aux.^2 .*(1 .-1 ./(1 .+exp.(-θ.*aux))))
end

function repeval(f:: Function)
    for i in 1:10000
        res = f(ones(n))
    end
end

@btime repeval(∇f)
@btime repeval(∇g)