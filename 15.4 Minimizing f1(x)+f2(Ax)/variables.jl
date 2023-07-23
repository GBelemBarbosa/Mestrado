using Distributions
using LinearAlgebra

#Dimensões & constantes
m=50
n=30
ϵ=eps() #Critério de parada
k_max=100 #Número máximo de iterações
ρ=1

#Variáveis aleatórias
A=randn(Float64, (m, n))
B=[-(i==j) for i=1:m, j=1:m]
c=zeros(m)

y₀=randn(Float64, m)