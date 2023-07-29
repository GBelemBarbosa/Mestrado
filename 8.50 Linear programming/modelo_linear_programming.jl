include("variables_linear_programming.jl")

f(x:: Vector{<:Number}; c=c)=c'*x

g(x:: Vector{<:Number}; A=A, b=b)=A*x.-b

oracle(λ:: Vector{<:Number}; A=A, b=b, c=c, n=n)=Float64.(I[1:n, argmin(c+A'*λ)])

γₖ(k:: Int64, gx:: Vector{<:Number})=1/sqrt(k+1) #Stepsize rule

include("../Métodos/Projected methods/dual_projected_subgradient_plot.jl")

x, p=dual_projected_subgradient(f, g, oracle, γₖ, λ₀, k_max, ϵ)