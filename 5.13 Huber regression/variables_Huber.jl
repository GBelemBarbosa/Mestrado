using Distributions
using LinearAlgebra
using Random

#Dimensões & constantes
m=80
n=100
out=5 #Number of outliers
σ=0.5 #STD comum => em Julia Normal(, σ)==N(, σ^2) em math notation
σₒ=5 #STD dos outliers
λ=1
α=λ #α-strongly convex
η=1
β=m*η+λ #β-smoothness (gradiente β-Lipschitz)
ϵ=eps() #Critério de parada caso tamanho do passo menor que ϵ
k_max=100 #Número máximo de iterações

#Variáveis aleatórias
xₒₚₜ=randn(n)
A=randn(Float64, (m, n))
y=A*xₒₚₜ
P=randperm(m) #Permutação para decidir quais serão outliers
y[P[1:out]].+=rand(Normal(0, σₒ), out) #Outliers
y[P[out+1:end]].+=rand(Normal(0, σ), m-out) 

x₀=zeros(n)