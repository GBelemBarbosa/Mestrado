using LinearAlgebra

include("variables.jl")

f₁(x:: Vector{<:Number})=norm(x, 2)^2/2

f₂(x:: Vector{<:Number})=norm(x, 1)

H(x:: Vector{<:Number}; A=A)=f₁(x)+f₂(A*x)

minimize_x(v:: Vector{<:Number}; A=A, ρ=ρ)=(ρ.*A'*A+I)\(ρ.*A'*v)

function proxf₁ρ(v:: Vector{<:Number}; f₁=f₁, ρ=ρ)
    aux=f₁(v)

    return (max(aux-ρ, 0)/aux).*v
end

Τ(λ:: Number, x:: Vector{<:Number})=[max(abs(x[i])-λ, 0)*sign(x[i]) for i=1:length(x)]

proxf₂ρ(v:: Vector{<:Number}; Τ=Τ, ρ=ρ)=Τ(ρ, v)

include("../Métodos/Lagrangian methods/alternating_direction_multipliers_plot.jl")

x, p₁=ADMM(H, A, B, c, minimize_x, proxf₂ρ, ρ, copy(y₀), k_max)

h₁(x:: Vector{<:Number})=0

h₂(z:: Vector{<:Number}, w:: Vector{<:Number})=f₁(w)+f₂(z)

A₂=vcat(A, I)
B₂=[-(i==j) for i=1:m+n, j=1:m+n]
c₂=zeros(m+n)

minimize_x₂(v:: Vector{<:Number}; A₂=A₂, ρ=ρ)=A₂'*A₂\(A₂'*v)

proxh₂ρ(v:: Vector{<:Number}; proxf₁ρ=proxf₁ρ, proxf₂ρ=proxf₂ρ, ρ=ρ, m=m)=vcat(proxf₂ρ(v[1:m]), proxf₁ρ(v[m+1:end]))

y₀₂=randn(Float64, m+n)

x, p₂=ADMM(H, A₂, B₂, c₂, minimize_x₂, proxh₂ρ, ρ, y₀₂, k_max)

α=ρ*opnorm(A, 2)^2 #ρ*λ_max(ATA)

β=ρ

proxf₁α(v:: Vector{<:Number}; α=α)=proxf₁ρ(v; ρ=α)

proxf₂β(v:: Vector{<:Number}; β=β)=Τ(β, v)

include("../Métodos/Lagrangian methods/alternating_direction_linearized_proximal_multipliers_plot.jl")

x, p₃=AD_LPMM(H, A, B, c, proxf₁α, proxf₂β, ρ, α, β, y₀, k_max)

plot(p₁, p₂, p₃)