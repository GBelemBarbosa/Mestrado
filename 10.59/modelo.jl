using LinearAlgebra
using Roots

include("variables.jl")

f(x:: Vector{<:Number})=0

h(x:: Vector{<:Number})=norm(x, 2)

function hμ(x:: Vector{<:Number}; μ=μ)
    nx=norm(x, 2)
    
    if nx>μ
        return nx-μ/2
    end

    return (nx^2)/(2*μ)
end

proxμh(x:: Vector{<:Number}; μ=μ)=(1-μ/max(norm(x, 2), μ)).*x

proxμhfalse(L:: Number, x:: Vector{<:Number}; proxμh=proxμh)=x.-proxμh(x) #Toma lugar de ∇Fμ para que y.-∇Fμy./L̃=proxμh(y) (com L̃=1 forçado)

PC(x:: Vector{<:Number})=min.(max.(x, -1), 1) #Projeção de x em C que é a caixa [-e, e]
    
include("../Métodos/Proximal methods/S-FISTA_plot.jl")

x, p=SFISTA(f, f, h, hμ, proxμhfalse, PC, 1-Lhμ, μ, α, x₀, k_max, eps()) #Lf=1-Lhμ=Lf-α/μ na chamada p/ q forçadamente L̃=Lf+α/μ=1 