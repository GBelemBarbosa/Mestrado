using Distributions
using LinearAlgebra
using StatsBase

#Dimensões & constantes
S=5 #Número de links
L=5 #Número de sources
ϵ=eps() #Critério de parada caso g(x) em norma infinito menor que ϵ
k_max=100 #Número máximo de iterações

#Variáveis aleatórias
M=abs.(randn(Float64, S)) #Máxima data rate para cada source 𝓈
cℓ=abs.(randn(Float64, L)) #Máxima capacidade para cada link ℓ
𝒮ℓ=[[sample(1:S, rand(1:S), replace=false)] for l=1:L] #Sources 𝓈 conectadas ao link ℓ
ℒ𝓈=[[l for l=1:L if s∈𝒮ℓ[l]] for s=1:S] #Links ℓ conectados a source 𝓈

λ₀=zeros(m)