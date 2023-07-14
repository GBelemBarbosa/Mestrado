using Distributions
using LinearAlgebra

#Dimensões & constantes
m=5
n=100
ϵ=eps() #Critério de parada caso g(x) em norma infinito menor que ϵ
k_max=100 #Número máximo de iterações

#Variáveis aleatórias
A=randn(Float64, (m, n))
b=randn(Float64, m)
c=randn(Float64, n)

λ₀=zeros(m)