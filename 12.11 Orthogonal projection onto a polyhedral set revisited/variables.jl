using Distributions
using LinearAlgebra

#Dimensões & constantes
p=30
n=50
ϵ=eps() #Critério de parada
k_max=20 #Número máximo de iterações
σ=1

#Variáveis aleatórias
A=randn(Float64, (p, n))
b=randn(Float64, p)
d=randn(Float64, n)

LF=p/σ

y₀=randn(Float64, p*n)