using Random
using StatsBase: sample
using LinearAlgebra
using Plots
using LaTeXStrings

include("../Group sparsity/group_sparse_functions.jl")

Random.seed!(1)

m       = 250
n       = 10*m
k_max   = 150
ϵ       = 10^-5
percent = 1

x⃰ = zeros(Float64, n) 
x⃰[sample(1:n, Int(m/10); replace=false)] .= 1
A  = randn(Float64, (m, n))
b  = A*x⃰.+0.01.*randn(Float64, m)

Lf = eigmax(A'A)
L  = 1.01*Lf
λ  = Lf/100
s  = round(Int64, percent*n)

f(x:: Array{<:Number}; A=A, b=b)  = norm(A*x.-b)^2/2
h(x:: Array{<:Number}; λ=λ)  = λ*norm(x, 0)
F(x:: Array{<:Number}; f=f, h=h)  = f(x)+h(x)
∇f(x:: Array{<:Number}; A=A, b=b) = A'*(A*x.-b)

Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = L, proxhL(L, x.-∂fx./L, s, λ, n)
pαₖ(αₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(αₖ, x.-∂fx./αₖ, s, λ, n)
proxα(α:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(1/α, x.-α.*∂fx, s, λ, n)
℘hλg(λₖ:: Number, k:: Int64, y:: Array{<:Number}, ∂fz:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(1/λₖ, y.-λₖ.*∂fz, s, λ, n)
Tλ(λₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(1/λₖ, x.-λₖ.*∂fx, s, λ, n)

x₀    = randn(n)
F⃰    = F(x⃰)
pltF  = plot([1, k_max], [F⃰, F⃰], label=L"F^\ast", title=L"l_0"*" penalty", xlabel="Iteration number", ylabel=L"F(x_k)", yscale=:log10)
pltnψ = plot(title=L"l_0"*" penalty", xlabel="Iteration number", ylabel=L"\|\|\psi(x_k)\|\|", yscale=:log10)

include("experiment.jl")