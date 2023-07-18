using Distributions
using LinearAlgebra

#Dimensões & constantes
n=30
p=15
ϵ=eps() #Critério de parada
k_max=100 #Número máximo de iterações

#Variáveis aleatórias
A=randn(Float64, (n, n))
A=A'*A
b=randn(Float64, n)

Lf₂=opnorm(A, 2) #Smooth constant de f no caso Euclidiano
Lf₁=maximum([maximum(abs.(A[:, i])) for i=1:n]) #Smooth constant de f no ℓ₁

x₀=randn(Float64, n)