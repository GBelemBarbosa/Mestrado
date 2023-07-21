using Distributions
using LinearAlgebra

#Dimensões & constantes
n=50
ϵ=eps() #Critério de parada
k_max=100 #Número máximo de iterações
λ=1 #Penalização
σ=1
LF=4/σ

#Variáveis aleatórias
d=reduce(vcat, rand(1:50).*ones(10) for i=1:n/10).+randn(Float64, n)./2

y₀=randn(Float64, n-1)