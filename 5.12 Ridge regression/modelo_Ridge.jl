include("variables_ridge.jl")

f(x:: Vector{Float64}; A=A, y=y, λ=λ)=(norm(A*x.-y)^2+λ*norm(x)^2)/2
#Outra alternativa que talvez seja mais estável numericamente e com tempos comparáveis em testes
#= 
function h(x:: Vector{Float64}; A=A, y=y, λ=λ)
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end 
=#

∇f(x:: Vector{Float64}; A=A, y=y, λ=λ)=A'*(A*x.-y).+λ.*x

tₖ(k:: Int64)=1/β #Stepsize rule p/ subgradient descent

include("../Métodos/subgradient_descent_plot.jl")

p₁=subgradient_descent(f, ∇f, tₖ, copy(x₀), k_max, ϵ)

include("../Métodos/acc_subgradient_descent_plot.jl")

p₂=acc_subgradient_descent(f, ∇f, β, copy(x₀), k_max, ϵ)

plot(p₁, p₂)