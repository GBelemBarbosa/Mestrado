using LinearAlgebra
using Plots
using LaTeXStrings
using KrylovKit
using Printf
using MAT
using SparseArrays

include("FISTA.jl")
include("BB.jl")
include("nmBB.jl")
include("proximal_subgradient.jl")

# Localização dos dados
lassodir = "../../Data-Lasso/"
problem  = "SC8"
vars     = matread(lassodir*problem*".mat")

const A, b, λ₁ = vars["A"], vec(vars["b"]), vars["lambda"] # λ₁ da norma l1 

const m, n    = size(A)
const k_max   = 300
const ϵ       = 10^-5 # Critério de parada para medida de estacionariedade ||ψ||

ATAx(x:: Array{<:Number}; A=A) = A'*(A*x)

const Lf = real(eigsolve(ATAx, n, 1, :LM, eltype(A))[1][1])
const L  = 1.01*Lf
const λ  = 0.1*norm(A'b, Inf)^2/(2*L) # Pode ser c*norm(A'b, Inf)^2/(2*L), onde 0<c<1
const x₀ = zeros(Float64, n)

f(x:: Array{<:Number}; A=A, b=b)    = norm(A*x.-b)^2/2
h(x:: Array{<:Number}; λ=λ)         = λ*norm(x, 0)
F(x:: Array{<:Number}; f=f, h=h)    = f(x)+h(x)
∇f(x:: Array{<:Number}; A=A, b=b)   = A'*(A*x.-b)

proxhL(L:: Number, x:: Vector{<:Number}; λ=λ) = [@inbounds x[i]*(abs(x[i])>sqrt(2*λ/L)) for i=eachindex(x)]

sλ    = @sprintf("%.2e", λ)
pltF  = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"F(x_k)", yscale=:log10, dpi=600)
pltnψ = plot(title=L"l_0"*" penalty, "*problem*", "*L"\lambda="*sλ, xlabel="CPU time", ylabel=L"\|\|\psi_k\|\|", yscale=:log10, dpi=600)

const α₀ = (sqrt(n)*10^-5)/norm(∇f(x₀).-∇f(x₀.+10^-5)) # Aproximação de 1/L entre x₀ e x₀.+10^-5
x_BB, histF, histnψ = BB(F, ∇f, proxhL, x₀, α₀, 10^-30, 10^30, k_max; ϵ=ϵ)   
plot!(pltF, histF, marker=:circle, label="SPG")
plot!(pltnψ, histnψ, marker=:circle, label="SPG")
println("nz(x_SPG)=", 1-norm(x_BB, 0)/n, ", f(x_SPG)=", f(x_BB))

# Não sei a razão, mas o primeiro nmBB rodado tem sempre um atraso, então rodo com 1 iteração e descarto
nmBB(F, ∇f, proxhL, x₀, 4, 0.01, 10^-30, 10^30, 1, 1; ϵ=ϵ)

x_mBB, histF, histnψ, ls_mBB = nmBB(F, ∇f, proxhL, x₀, α₀, 4, 0.01, 10^-30, 10^30, 1, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:circle, label="SPGmLS")
plot!(pltnψ, histnψ, marker=:circle, label="SPGmLS")
println("nz(x_SPGmLS)=", 1-norm(x_mBB, 0)/n, ", f(x_SPGmLS)=", f(x_mBB), ", #LS=", ls_mBB)

x_nmBB, histF, histnψ, ls_nmBB = nmBB(F, ∇f, proxhL, x₀, α₀, 4, 0.01, 10^-30, 10^30, 5, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:circle, label="SPGnmLS")
plot!(pltnψ, histnψ, marker=:circle, label="SPGnmLS")
println("nz(x_SPGnmLS)=", 1-norm(x_nmBB, 0)/n, ", f(x_SPGnmLS)=", f(x_nmBB), ", #LS=", ls_nmBB)

x_PG, histF, histnψ = PG(F, ∇f, proxhL, x₀, L, k_max; ϵ=ϵ)
plot!(pltF, histF, linestyle=:dash, label="PG")
plot!(pltnψ, histnψ, linestyle=:dash, label="PG")
println("nz(x_PG)=", 1-norm(x_PG, 0)/n, ", f(x_PG)=", f(x_PG))

x_FISTA, histF, histnψ = FISTA(F, ∇f, proxhL, x₀, L, k_max; ϵ=ϵ)
plot!(pltF, histF, linestyle=:dash, label="FISTA")
plot!(pltnψ, histnψ, linestyle=:dash, label="FISTA")
println("nz(x_FISTA)=", 1-norm(x_FISTA, 0)/n, ", f(x_FISTA)=", f(x_FISTA))

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_part.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_part.png")

#savefig(pltF, "Experiments/Plots/F"*problem*"_"*sλ*"_full.png")
#savefig(pltnψ, "Experiments/Plots/nPhi"*problem*"_"*sλ*"_full.png")

plot(pltF, pltnψ)