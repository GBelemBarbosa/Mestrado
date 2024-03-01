using Random
using StatsBase: sample
using LinearAlgebra
using Plots
using LaTeXStrings

include("../Group sparsity/group_sparse_functions.jl")

Random.seed!(1)

m       = 250
n       = 10*m
k_max   = 60
ϵ       = 10^-5

x⃰ = zeros(Float64, n) 
x⃰[sample(1:n, Int(m/10); replace=false)] .= 1
A  = randn(Float64, (m, n))
b  = A*x⃰.+0.01.*randn(Float64, m)

Lf = eigmax(A'A)
L  = 1.01*Lf
λ  = m/1000

f(x:: Array{<:Number}; A=A, b=b)  = norm(A*x.-b)^2/2
h(x:: Array{<:Number}; λ=λ, p=p)  = λ*norm(x, 1)
F(x:: Array{<:Number}; f=f, h=h)  = f(x)+h(x)
∇f(x:: Array{<:Number}; A=A, b=b) = A'*(A*x.-b)

Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL1=proxhL1, n=n, λ=λ) = L, proxhL1(L, x.-∂fx./L, λ)
pαₖ(αₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL1=proxhL1, n=n, λ=λ) = proxhL1(αₖ, x.-∂fx./αₖ, λ)
proxα(α:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL1=proxhL1, n=n, λ=λ) = proxhL1(1/α, x.-α.*∂fx, λ)
℘hλg(αₖ:: Number, λₖ:: Number, k:: Int64, y:: Array{<:Number}, ∂fz:: Array{<:Number}; proxhL1=proxhL1, n=n, λ=λ, s=s) = proxhL1(1/λₖ, y.-λₖ.*∂fz, λ)
Tλ(λₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL1=proxhL1, n=n, λ=λ) = proxhL1(1/λₖ, x.-λₖ.*∂fx, λ)

x₀  = randn(n)
F⃰  = F(x⃰)
pltF  = plot([1, k_max], [F⃰, F⃰], label=L"F_{opt}", title=L"l_1"*" penalty", xlabel="Iteration number", ylabel=L"F(x_k)", yscale=:log10)
pltnψ = plot(title=L"l_0"*" penalty", xlabel="Iteration number", ylabel=L"\|\|\psi(x_k)\|\|", yscale=:log10)

include("experiment.jl")