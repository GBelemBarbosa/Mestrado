using LinearAlgebra

include("variables.jl")

f(x:: Vector{<:Number}; d=d)=norm(x.-d, 2)^2/2

g(x:: Vector{<:Number})=0

step(x:: Vector{<:Number}; d=d)=x.+d

Lₖ(L:: Number, k:: Int64, y:: Vector{<:Number}, Ax:: Vector{<:Number}; b=b)=L, min.(Ax.-L*y, b)

include("../Métodos/Proximal methods/dual_proximal_subgradient_plot.jl")

p₁=dual_proximal_subgradient(f, g, step, A, Lₖ, copy(y₀), LF, k_max, ϵ)

include("../Métodos/Proximal methods/fast_dual_proximal_subgradient_plot.jl")

p₂=fast_dual_proximal_subgradient(f, g, step, A, Lₖ, y₀, LF, k_max, ϵ)

plot(p₁, p₂)