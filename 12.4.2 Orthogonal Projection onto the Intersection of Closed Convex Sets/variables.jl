using Distributions
using LinearAlgebra

#Dimensões & constantes
n=10
p=3
ϵ=eps() #Critério de parada
k_max=20 #Número máximo de iterações
σ=1
LF=p/σ

#Variáveis aleatórias
d=randn(Float64, n)

y₀=randn(Float64, n*p)