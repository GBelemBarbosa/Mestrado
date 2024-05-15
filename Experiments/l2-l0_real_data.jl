using LinearAlgebra
using Plots
using LaTeXStrings
using KrylovKit
using Printf
using MAT
using SparseArrays

include("../Group sparsity/group_sparse_functions.jl")
include("../Métodos/Proximal methods hist/FISTA.jl")

# Localização dos dados e métodos
lassodir = "../../Data-Lasso/"
problem  = "SC8"
vars     = matread(lassodir*problem*".mat")

methods = ["SPG", "SPGnmLS", "SPGmLS", "PG", "FISTA"]

const A, b, λ1 = vars["A"], vec(vars["b"]), vars["lambda"] # λ1 da norma l1 

const m, n    = size(A)
const k_max   = 300
const ϵ       = 10^-5
const percent = 1

ATAx(x:: Array{<:Number}; A=A) = A'*(A*x)

const Lf = real(eigsolve(ATAx, n, 1, :LM, eltype(A))[1][1])
const L  = 1.01*Lf
const λ  = 0.1*norm(A'b, Inf)^2/(2*L) # Pode ser c*norm(A'b, Inf)^2/(2*L), onde 0<c<1
const x₀ = zeros(Float64, n)

f(x:: Array{<:Number}; A=A, b=b)    = norm(A*x.-b)^2/2
h(x:: Array{<:Number}; λ=λ)         = λ*norm(x, 0)
h1(x:: Array{<:Number}; λ1=λ1)      = λ1*norm(x, 1)
F(x:: Array{<:Number}; f=f, h=h)    = f(x)+h(x)
F1(x:: Array{<:Number}; f=f, h1=h1) = f(x)+h1(x)
∇f(x:: Array{<:Number}; A=A, b=b)   = A'*(A*x.-b)

Lₖl1(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ1=λ1) = L, proxhL1(L, x.-∂fx./L, λ1)

x_FISTA, histF, histnψ = FISTA(F1, ∇f, Lₖl1, x₀, L, k_max; ϵ=ϵ)
pltFl1  = plot(histF, linestyle=:dash, label="FISTA", title=L"l_1"*" penalty, "*problem, xlabel="CPU time", ylabel=L"F(x_k)", yscale=:log10)
pltnψl1 = plot(histnψ, linestyle=:dash, label="FISTA", title=L"l_1"*" penalty, "*problem, xlabel="CPU time", ylabel=L"\|\|\psi(x_k)\|\|", yscale=:log10)
println("l1 FISTA: ", 1-norm(x_FISTA, 0)/n, " & ", f(x_FISTA))

Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = L, proxhL(L, x.-∂fx./L, λ)
pαₖ(αₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(αₖ, x.-∂fx./αₖ, λ)
proxα(α:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/α, x.-α.*∂fx, λ)
℘hλg(λₖ:: Number, k:: Int64, y:: Array{<:Number}, ∂fz:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/λₖ, y.-λₖ.*∂fz, λ)
Tλ(λₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/λₖ, x.-λₖ.*∂fx, λ)

sλ    = @sprintf("%.2e", λ)
pltF  = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"F(x_k)", yscale=:log10, dpi=600)
pltnψ = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"\|\|\psi_k\|\|", yscale=:log10, dpi=600)

include("experiment.jl")

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_part.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_part.png")

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_full.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_full.png")