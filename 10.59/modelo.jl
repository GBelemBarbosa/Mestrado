using LinearAlgebra
using Roots

include("variables.jl")

f(x:: Vector{Float64})=0

h(x:: Vector{Float64})=norm(x, 2)

function hμ(x:: Vector{Float64}; μ=μ)
    nx=norm(x, 2)
    
    if nx>μ
        return nx-μ/2
    end

    return (nx^2)/(2*μ)
end

proxμh(x:: Vector{Float64}; μ=μ)=(1-μ/max(norm(x, 2), μ)).*x

proxμhfalse(x:: Vector{Float64}; L=L, proxμh=proxμh)=x.-proxμh(x) #Toma lugar de ∇f para que x.-∇fx./L=proxμh com L=1 (incorreto) na chamada

PC(x:: Vector{Float64})=min.(max.(x, -1), 1) #Projeção de x em C que é a caixa [-e, e]
    
include("../Métodos/Proximal methods/S-FISTA_plot.jl")

p=SFISTA(f, f, h, hμ, proxμhfalse, PC, 1, μ, α, x₀, k_max, eps())