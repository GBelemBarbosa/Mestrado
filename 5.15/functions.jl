include("variables.jl")

f(x:: Vector{Float64}; A=A, y=y, λ=λ)=(norm(A*x-y)^2+λ*norm(x)^2)/2
#Outra alternativa que talvez seja mais estável numericamente e com tempos comparáveis em testes
#= function h(x:: Vector{Float64}; A=A, y=y, λ=λ)
    aux=A*x-y
    return (aux⋅aux+λ*x⋅x)/2
end =#

∇f(x:: Vector{Float64}; A=A, y=y, λ=λ)=A'*(A*x.-y).+λ.*x

include("gradient_descent_plot.jl")

gradient_descent(f, ∇f, β, x₀, t_max)