using Distributions
using LinearAlgebra

#Dimensões & constantes
m=5
n=10
ϵ=eps() #Critério de parada caso tamanho do passo (em termos de distância à S₁) menor que ϵ
k_max=100 #Número máximo de iterações

#Variáveis aleatórias
xₒₚₜ=randn(n)
A=randn(Float64, (m, n))
b=randn(Float64, m)

x₀=zeros(n)