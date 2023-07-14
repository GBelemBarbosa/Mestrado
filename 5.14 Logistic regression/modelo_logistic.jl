include("variables_logistic.jl")

f(θ:: Vector{Float64}; X=X, y=y, λ=λ)=sum(log.(1 .+exp.(-y.*X*θ)))+λ*norm(θ)^2/2

function ∇f(θ:: Vector{Float64}; X=X, y=y, λ=λ)
    aux=X'y
    return λ.*θ.-aux./(1 .+exp.(aux.*θ))
end

tₖ(k:: Int64, ∂f:: Vector{Float64})=1/β #Stepsize rule p/ subgradient descent

include("../Métodos/Descent methods/subgradient_descent_plot.jl")

p₁=subgradient_descent(f, ∇f, tₖ, copy(θ₀), k_max, ϵ)

include("../Métodos/Descent methods/acc_subgradient_descent_plot.jl")

p₂=acc_subgradient_descent(f, ∇f, β, copy(θ₀), k_max, ϵ)

plot(p₁, p₂)