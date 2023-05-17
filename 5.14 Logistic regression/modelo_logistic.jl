include("variables_logistic.jl")

f(θ:: Vector{Float64}; X=X, y=y, λ=λ)=sum(log.(1 .+exp.(-y.*X*θ)))+λ*norm(θ)^2/2

function ∇f(θ:: Vector{Float64}; X=X, y=y, λ=λ)
    aux=X'y
    return λ.*θ.-aux./(1 .+exp.(aux.*θ))
end

include("../Métodos/gradient_descent_plot.jl")

p₁=gradient_descent(f, ∇f, β, θ₀, t_max)

include("../Métodos/acc_gradient_descent_plot.jl")

p₂=acc_gradient_descent(f, ∇f, β, θ₀, t_max)

plot(p₁, p₂)