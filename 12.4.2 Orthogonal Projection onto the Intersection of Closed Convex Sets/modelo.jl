using LinearAlgebra

include("variables.jl")

f(x:: Vector{Float64}; d=d)=norm(x.-d, 2)^2/2

g(x:: Vector{Float64})=0

step(x:: Vector{Float64}; d=d)=x.+d

Α(x:: Vector{Float64}; p=p)=reduce(vcat, x for i=1:p)

ΑT(x:: Vector{Float64}; p=p)=sum(x[i*n+1:(i+1)*n] for i=0:p-1)

PC(x:: Vector{Float64}; n=n)=[min(x[1], 0); x[2:n]; x[n+1:2*n-1]; max(x[2*n], 0); x[2*n+1:end]./norm(x[2*n+1:end], 2)] #=Projeção no set 
cuja primeira coordenada é não negativa, cuja última é não negativa, e na esfera unitária, respectivamente=#

Lₖ(L:: Number, k:: Int64, y:: Vector{Float64}, Αx:: Vector{Float64}; b=b, PC=PC)=L, PC(Αx.-L*y)

include("../Métodos/Proximal methods/dual_proximal_subgradient_plot.jl")

p₁=dual_proximal_subgradient(f, g, step, Α, ΑT, Lₖ, y₀, LF, k_max, ϵ)

include("../Métodos/Proximal methods/fast_dual_proximal_subgradient_plot.jl")

p₂=fast_dual_proximal_subgradient(f, g, step, Α, ΑT, Lₖ, y₀, LF, k_max, ϵ)

plot(p₁, p₂)