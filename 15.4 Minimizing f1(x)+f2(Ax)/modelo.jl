using LinearAlgebra

include("variables.jl")

f₁(x:: Vector{Float64})=norm(x, 2)^2/2

f₂(x:: Vector{Float64})=norm(x, 1)

H(x:: Vector{Float64}; A=A)=f₁(x)+f₂(A*x)

minimize_x(v:: Vector{Float64}; A=A, ρ=ρ)=(ρ.*A'*A+I)\(ρ.*A'*v)

function proxf₁ρ(v:: Vector{Float64}; f₁=f₁, ρ=ρ)
    aux=f₁(v)

    return (max(aux-ρ, 0)/aux).*v
end

Τ(λ:: Number, x:: Vector{Float64})=[max(abs(x[i])-λ, 0)*sign(x[i]) for i=1:length(x)]

proxf₂ρ(v:: Vector{Float64}; Τ=Τ, ρ=ρ)=Τ(ρ, v)

include("../Métodos/Lagrangian methods/alternating_direction_multipliers_plot.jl")

p₁=ADMM(H, A, B, c, minimize_x, proxf₂ρ, ρ, copy(y₀), k_max, ϵ)

h₁(x:: Vector{Float64})=0

h₂(z:: Vector{Float64}, w:: Vector{Float64})=f₁(w)+f₂(z)

A₂=vcat(A, I)
B₂=[-(i==j) for i=1:m+n, j=1:m+n]
c₂=zeros(m+n)

minimize_x₂(v:: Vector{Float64}; A₂=A₂, ρ=ρ)=A₂'*A₂\(A₂'*v)

proxh₂ρ(v:: Vector{Float64}; proxf₁ρ=proxf₁ρ, proxf₂ρ=proxf₂ρ, ρ=ρ, m=m)=vcat(proxf₂ρ(v[1:m]), proxf₁ρ(v[m+1:end]))

y₀₂=randn(Float64, m+n)

p₂=ADMM(H, A₂, B₂, c₂, minimize_x₂, proxh₂ρ, ρ, y₀₂, k_max, ϵ)

α=ρ*opnorm(A, 2)^2 #ρ*λmax(ATA)

β=ρ

proxf₁α(v:: Vector{Float64}; α=α)=proxf₁ρ(v; ρ=α)

proxf₂β(v:: Vector{Float64}; β=β)=Τ(β, v)

include("../Métodos/Lagrangian methods/alternating_direction_linearized_proximal_multipliers_plot.jl")

p₃=AD_LPMM(H, A, B, c, proxf₁α, proxf₂β, ρ, α, β, y₀, k_max, ϵ)

plot(p₁, p₂, p₃)