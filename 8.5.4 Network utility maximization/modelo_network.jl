include("variables_network.jl")

r=randn(Float64, S) #Coeficientes que definem utility rate uₛ (linear e portanto concava) de cada source s

f(x:: Vector{Float64}; c=c)=-r'*x

g(x:: Vector{Float64}; cℓ=cℓ, 𝒮ℓ=𝒮ℓ, L=L)=[sum(x[i for i=𝒮ℓ(l)])-cℓ[l] for l=1:L]

oracle(λ:: Vector{Float64}; r=r, M=M, ℒ𝓈=ℒ𝓈, S=S)=[max(r(s)*M(s), 0) for s=1:S] 

γₖ(k:: Int64, ∂f:: Vector{Float64})=1/sqrt(k+1) #Stepsize rule

include("../Métodos/Projected methods/dual_projected_subgradient_plot.jl")

#p=dual_projected_subgradient(f, g, oracle, γₖ, λ₀, k_max, ϵ)