using Distributions
using LinearAlgebra

#Dimensões & constantes
ϵ=10^(-2) #Critério de tolerância
k_max=100 #Número máximo de iterações
ℓₕ=1 #Lipshitz constant de h
μ=ϵ/ℓₕ^2 #Smoothing parameter p/ obter convergência em 𝛰(1/ϵ) iterações
α=1 #Smooth approximation parameter
Lhμ=α/μ #Smooth constant de hμ

x₀=randn(Float64, 10)