using LinearAlgebra

include("variables.jl")

f(x:: Vector{<:Number}; d=d)=norm(x.-d, 2)^2/2

g(x:: Vector{<:Number})=0

step(x:: Vector{<:Number}; d=d)=x.+d

Α(x:: Vector{<:Number}; p=p)=reduce(vcat, x for i=1:p)

ΑT(x:: Vector{<:Number}; p=p)=sum(x[i*n+1:(i+1)*n] for i=0:p-1)

normal_ai=reduce(vcat, A[i, :]./norm(A[i, :]) for i=1:p)

PC(x:: Vector{<:Number}; A=A, b=b, normal_ai=normal_ai, n=n)=reduce(vcat, max(A[i, :]'*x[(i-1)*n+1:i*n]-b[i], 0).*normal_ai[(i-1)*n+1:i*n]  for i=1:p) 

Lₖ(L:: Number, k:: Int64, y:: Vector{<:Number}, Αx:: Vector{<:Number}; b=b, PC=PC)=L, PC(Αx.-L*y)

include("../Métodos/Proximal methods/dual_proximal_subgradient_plot.jl")

x, p₁=dual_proximal_subgradient(f, g, step, Α, ΑT, Lₖ, copy(y₀), LF, k_max, ϵ)

include("../Métodos/Proximal methods/fast_dual_proximal_subgradient_plot.jl")

x, p₂=fast_dual_proximal_subgradient(f, g, step, Α, ΑT, Lₖ, y₀, LF, k_max, ϵ)

plot(p₁, p₂)