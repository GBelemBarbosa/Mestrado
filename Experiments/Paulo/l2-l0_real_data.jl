using LinearAlgebra
using Plots
using LaTeXStrings
using KrylovKit
using Printf
using MAT
using SparseArrays

include("group_sparse_functions.jl")
include("FISTA.jl")
include("BB.jl")
include("nmBB.jl")

# Localização dos dados
lassodir = "../../Data-Lasso/"
problem  = "SC12"
vars     = matread(lassodir*problem*".mat")

const A, b, λ₁ = vars["A"], vec(vars["b"]), vars["lambda"] # λ₁ da norma l1 

const m, n    = size(A)
const k_max   = 300
const ϵ       = 10^-5 # Critério de parada para medida de estacionariedade ||ψ||

ATAx(x:: Array{<:Number}; A=A) = A'*(A*x)

const Lf = real(eigsolve(ATAx, n, 1, :LM, eltype(A))[1][1])
const L  = 1.01*Lf
const λ  = 0.1*norm(A'b, Inf)/L # Pode ser c*norm(A'b, Inf)/L, onde 0<c<1
const x₀ = zeros(Float64, n)

f(x:: Array{<:Number}; A=A, b=b)    = norm(A*x.-b)^2/2
h(x:: Array{<:Number}; λ=λ)         = λ*norm(x, 0)
F(x:: Array{<:Number}; f=f, h=h)    = f(x)+h(x)
∇f(x:: Array{<:Number}; A=A, b=b)   = A'*(A*x.-b)

proxhL(L:: Number, x:: Array{<:Number}; proxl0L=proxl0L, λ=λ) = proxl0L(L, x, λ)

sλ    = @sprintf("%.2e", λ)
pltF  = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"F(x_k)", yscale=:log10, dpi=600)
pltnψ = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"\|\|\psi(x_k)\|\|", yscale=:log10, dpi=600)

const α₀ = norm(∇f(x₀).-∇f(x₀.+10^-5))/(sqrt(n)*10^-5) # Aproximação de L entre x₀ e x₀.+10^-5
x_BB, histF, histnψ = BB(F, ∇f, proxhL, x₀, α₀, 10^-30, 10^30, 5, k_max; ϵ=ϵ)   
plot!(pltF, histF, marker=:circle, label="BB")
plot!(pltnψ, histnψ, marker=:circle, label="BB")
println("nz(x_BB)=", 1-norm(x_BB, 0)/n, ", f(x_BB)=", f(x_BB))

x_mBB, histF, histnψ, ls_mBB = nmBB(F, ∇f, proxhL, x₀, 4, 0.01, 10^-30, 10^30, 1, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:circle, label="mBB")
plot!(pltnψ, histnψ, marker=:circle, label="mBB")
println("nz(x_mBB)=", 1-norm(x_mBB, 0)/n, ", f(x_mBB)=", f(x_mBB), ", #LS=", ls_mBB)

x_nmBB, histF, histnψ, ls_nmBB = nmBB(F, ∇f, proxhL, x₀, 4, 0.01, 10^-30, 10^30, 5, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:circle, label="nmBB")
plot!(pltnψ, histnψ, marker=:circle, label="nmBB")
println("nz(x_nmBB)=", 1-norm(x_nmBB, 0)/n, ", f(x_nmBB)=", f(x_nmBB), ", #LS=", ls_nmBB)

x_FISTA, histF, histnψ = FISTA(F, ∇f, proxhL, x₀, L, k_max; ϵ=ϵ)
plot!(pltF, histF, linestyle=:dash, label="FISTA")
plot!(pltnψ, histnψ, linestyle=:dash, label="FISTA")
println("nz(x_FISTA)=", 1-norm(x_FISTA, 0)/n, ", f(x_FISTA)=", f(x_FISTA))

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_part.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_part.png")

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_full.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_full.png")