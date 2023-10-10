include("variables_logistic.jl")

f(θ:: Vector{<:Number}; X=X, y=y, λ=λ)=sum(log.(1 .+exp.(-y.*X*θ)))+λ*norm(θ)^2/2

function ∇f(θ:: Vector{<:Number}; X=X, y=y, λ=λ)
    aux=X'y
    return λ.*θ.-aux./(1 .+exp.(aux.*θ))
end

tₖ(k:: Int64, ∂f:: Vector{<:Number})=1/β #Stepsize rule p/ subgradient descent

include("../Métodos/Descent methods/subgradient_descent_plot.jl")

θ, p₁=subgradient_descent(f, ∇f, tₖ, copy(θ₀), k_max)

include("../Métodos/Descent methods/acc_subgradient_descent_plot.jl")

θ, p₂=acc_subgradient_descent(f, ∇f, β, θ₀, k_max)

plot(p₁, p₂)