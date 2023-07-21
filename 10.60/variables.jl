using Distributions
using LinearAlgebra

#Dimensões & constantes
m=20
n=30
p=15
ϵ=10^(-2) #Critério de tolerância
k_max=100 #Número máximo de iterações
λ=abs(randn(Float64)) #Penalização

#Variáveis aleatórias
A=randn(Float64, (m, n))
b=randn(Float64, m)
D=randn(Float64, (p, n))

Lf=opnorm(A, 2)^2 #Smooth constant de f
nD=opnorm(D, 2)
sp=sqrt(p)
μ=2*nD*ϵ/(nD*sp+sqrt(p*nD^2+2*ϵ*Lf))/sp #Smoothing parameter p/ obter convergência em 𝛰(1/ϵ) iterações
α=1 #Smooth approximation parameter
Lhμ=α/μ #Smooth constant de hμ
L̃=Lf+Lhμ #Smooth constant de Fμ

x₀=randn(Float64, n)