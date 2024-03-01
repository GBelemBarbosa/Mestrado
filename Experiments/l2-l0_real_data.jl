using LinearAlgebra
using Plots
using LaTeXStrings
using KrylovKit
using Printf
using MAT
using SparseArrays

include("../Group sparsity/group_sparse_functions.jl")
include("../Métodos/Proximal methods hist/FISTA.jl")

# Localização dos dados
lassodir = "../../Data-Lasso/"
problem  = "SC12"
vars     = matread(lassodir*problem*".mat")

const A, b, λ1 = vars["A"], vec(vars["b"]), vars["lambda"] # λ1 da norma l1 

const m, n    = size(A)
const k_max   = 300
const ϵ       = 10^-5
const percent = 1

ATAx(x:: Array{<:Number}; A=A) = A'*(A*x)

const Lf = real(eigsolve(ATAx, n, 1, :LM, eltype(A))[1][1])
const L  = 1.01*Lf
const λ  = 0.1*norm(A'b, Inf)/L # Pode ser c*norm(A'b, Inf)/L, onde 0<c<1
const s  = round(Int64, percent*n)
const x₀ = zeros(Float64, n)

f(x:: Array{<:Number}; A=A, b=b)    = norm(A*x.-b)^2/2
h(x:: Array{<:Number}; λ=λ)         = λ*norm(x, 0)
h1(x:: Array{<:Number}; λ1=λ1)      = λ1*norm(x, 1)
F(x:: Array{<:Number}; f=f, h=h)    = f(x)+h(x)
F1(x:: Array{<:Number}; f=f, h1=h1) = f(x)+h1(x)
∇f(x:: Array{<:Number}; A=A, b=b)   = A'*(A*x.-b)

Lₖl1(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ1=λ1, s=s) = L, proxhL1(L, x.-∂fx./L, λ1)

x_FISTA, histF, histnψ = FISTA(F1, ∇f, Lₖl1, x₀, L, k_max; ϵ=ϵ)
pltFl1  = plot(histF, linestyle=:dash, label="FISTA", title=L"l_1"*" penalty, "*problem, xlabel="CPU time", ylabel=L"F(x_k)", yscale=:log10)
pltnψl1 = plot(histnψ, linestyle=:dash, label="FISTA", title=L"l_1"*" penalty, "*problem, xlabel="CPU time", ylabel=L"\|\|\psi(x_k)\|\|", yscale=:log10)
println("l1 FISTA: ", 1-norm(x_FISTA, 0)/n, " & ", f(x_FISTA))

Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = L, proxhL(L, x.-∂fx./L, s, λ, n)
pαₖ(αₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(αₖ, x.-∂fx./αₖ, s, λ, n)
proxα(α:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(1/α, x.-α.*∂fx, s, λ, n)
℘hλg(λₖ:: Number, k:: Int64, y:: Array{<:Number}, ∂fz:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(1/λₖ, y.-λₖ.*∂fz, s, λ, n)
Tλ(λₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ=λ, s=s) = proxhL(1/λₖ, x.-λₖ.*∂fx, s, λ, n)

sλ    = @sprintf("%.2e", λ)
pltF  = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"F(x_k)", yscale=:log10, dpi=600)
pltnψ = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"\|\|\psi(x_k)\|\|", yscale=:log10, dpi=600)

include("experiment.jl")

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_part.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_part.png")

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_full.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_full.png")