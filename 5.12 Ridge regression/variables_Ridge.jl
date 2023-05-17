using Distributions
using LinearAlgebra

#Dimensões & constantes
m=80
n=100
σ=0.5 #STD => em Julia Normal(, σ)==N(, σ^2) em math notation
λ=1
α=λ #α-strongly convex
β=4(m+n)+λ #β-smoothness (gradiente β-Lipschitz)
ϵ=eps() #Critério de parada caso tamanho do passo menor que ϵ
t_max=50 #Número máximo de iterações

#Variáveis aleatórias
xₒₚₜ=randn(n)
A=randn(Float64, (m, n))
y=A*xₒₚₜ.+rand(Normal(0, σ), m)

x₀=zeros(n)