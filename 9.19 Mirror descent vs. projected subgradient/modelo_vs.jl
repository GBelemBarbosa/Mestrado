using LinearAlgebra
using Roots

include("variables_vs.jl")

f(x:: Vector{Float64}; A=A, b=b)=norm(A*x.-b, 1)

∂f(x:: Vector{Float64}; A=A, b=b)=A'*sign.(A*x.-b)

aux(x:: Vector{Float64}, μ:: T) where T=sum(max.(x.-μ, 0))-1

PC(x:: Vector{Float64}; aux=aux)=max.(x.-find_zero((μ -> aux(x, μ)), (maximum(x), min(minimum(x), 0)-1)), 0) #Projeção de x em C
    
tₖ(k:: Int64, ∂f:: Vector{Float64})=sqrt(2/(k+1))/norm(∂f, 2) #Stepsize rule p/ projected subgradient

include("../Métodos/Projected methods/projected_subgradient_plot.jl")

p₁=projected_subgradient(f, ∂f, PC, tₖ, x₀, k_max, ϵ)

∇ω(x:: Vector{Float64})=0

function step(x:: Vector{Float64}, ∂:: Vector{Float64}; n=n)
    divergence=0.0
    x_new=Float64[]
    
    for i=1:n
        aux=x[i]*exp(-∂[i])
        push!(x_new, aux)
        divergence+=aux
    end

    return x_new./divergence
end

dual_norm(x:: Vector{Float64})=norm(x, Inf)

tₖ(k:: Int64, ∂f:: Vector{Float64})=sqrt(2/(k+1))/norm(∂f, Inf) #Stepsize rule p/ mirror descent

include("../Métodos/Non-Euclidean Methods/mirror_descent_plot.jl")

p₂=mirror_descent(f, ∂f, ∇ω, step, tₖ, dual_norm, x₀, k_max, ϵ)

plot(p₁, p₂)