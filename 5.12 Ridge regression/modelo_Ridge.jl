include("variables_ridge.jl")

f(x:: Vector{<:Number}; A=A, y=y, λ=λ)=(norm(A*x.-y)^2+λ*norm(x)^2)/2
#Outra alternativa que talvez seja mais estável numericamente e com tempos comparáveis em testes
#= 
function h(x:: Vector{<:Number}; A=A, y=y, λ=λ)
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end 
=#

∇f(x:: Vector{<:Number}; A=A, y=y, λ=λ)=A'*(A*x.-y).+λ.*x

tₖ(k:: Int64, ∂f:: Vector{<:Number})=1/β #Stepsize rule p/ subgradient descent

include("../Métodos/Descent methods/subgradient_descent_plot.jl")

x, p₁=subgradient_descent(f, ∇f, tₖ, copy(x₀), k_max)

include("../Métodos/Descent methods/acc_subgradient_descent_plot.jl")

x, p₂=acc_subgradient_descent(f, ∇f, β, x₀, k_max)

plot(p₁, p₂)