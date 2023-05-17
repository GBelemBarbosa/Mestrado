using Distributions
using LinearAlgebra
using Random

#Dimensões & constantes
m=100
n=50
λ=1
α=λ #α-strongly convex
β=100 #β-smoothness (gradiente β-Lipschitz)
ϵ=eps() #Critério de parada caso tamanho do passo menor que ϵ
t_max=100 #Número máximo de iterações

#Variáveis aleatórias
θₒₚₜ=ones(n)
X=randn(Float64, (m, n))
z=X*θₒₚₜ
p=1 ./(1 .+exp.(-z)) #Vetor de probabilidades
y=2 .*rand.(Bernoulli.(p)).-1
P=randperm(m) #Permutação para decidir quais serão outliers

θ₀=randn(n)