using LinearAlgebra

include("variables.jl")

f(x:: Vector{<:Number}; d=d)=norm(x.-d, 2)^2/2

g(x:: Vector{<:Number}; λ=λ)=λ*norm(x, 1)

step(x:: Vector{<:Number}; d=d)=x.+d

A=[(i==j)-(i+1==j) for i=1:n-1, j=1:n]

Τ(λ:: Number, x:: Vector{<:Number})=[max(abs(x[i])-λ, 0)*sign(x[i]) for i=1:length(x)] #Soft threshholding

Lₖ(L:: Number, k:: Int64, y:: Vector{<:Number}, Ax:: Vector{<:Number}; b=b, Τ=Τ)=L, Τ(L*λ, Ax.-L*y)

include("dual_proximal_subgradient_plot.jl")

x⃰, p₁=dual_proximal_subgradient(f, g, step, A, Lₖ, y₀, LF, k_max, ϵ)

include("fast_dual_proximal_subgradient_plot.jl")

x⃰f, p₂=fast_dual_proximal_subgradient(f, g, step, A, Lₖ, y₀, LF, k_max, ϵ)

plot(p₁, p₂)

p=scatter(eachindex(d), d, label=L"\textbf{d}")

p₃=scatter(p, eachindex(x⃰), x⃰, label=L"\textbf{x*}")

p₄=scatter(p, eachindex(x⃰f), x⃰f, label=L"\textbf{x*}\ (fast)")

plot(p₃, p₄)