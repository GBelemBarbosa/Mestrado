using Distributions
using LinearAlgebra

#Dimensões & constantes
m=80
n=100
σ=0.25
λ=1.0
α=λ
β=4(m+n)+λ
ϵ=10^-10
t_max=50

#Variáveis aleatórias
xₒₚₜ=randn(n)
A=randn(Float64, (m, n))
y=A*xₒₚₜ.+rand(Normal(0, σ), m)

x₀=rand(n)