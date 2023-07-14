using Distributions
using LinearAlgebra

#Dimensões & constantes
m=50
n=50
ϵ=eps() #Critério de parada caso g(x) em norma infinito menor que ϵ
k_max=100 #Número máximo de iterações
Θ=1 #limite superior na metade do quadrado do diâmetro de C, que é o unit simplex

#Variáveis aleatórias
A=randn(Float64, (m, n))
b=randn(Float64, m)

x₀=randn(Float64, n)