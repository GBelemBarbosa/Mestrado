using LinearAlgebra
using Plots
using LaTeXStrings
using KrylovKit
using Printf
using MAT
using SparseArrays
using BenchmarkProfiles
using FileIO, JLD2

include("../Group sparsity/group_sparse_functions.jl")
include("../Métodos/Proximal methods hist/FISTA.jl")
include("experiment_performance.jl")

methods = ["SPG", "SPGmLS", "SPGnmLS", "mAPG", "nmAPG", "PG", "FISTA", "newAPG_vs", "newAPG_vs_x_2", "newAPG_vs_y_2"]

const k_max   = 5000
const ϵ₀      = 10^-5

converg_T = Array{Float64}(undef, 25, length(methods))
F_hist    = Array{Float64}(undef, 25, length(methods))

for i=1:25
    # Localização dos dados
    lassodir = "../../Data-Lasso/"
    problem  = "SC"*string(i)
    vars     = matread(lassodir*problem*".mat")

    println(problem*":")

    A, b, λ1 = vars["A"], vec(vars["b"]), vars["lambda"] # λ1 da norma l1 
    m, n     = size(A)

    ATAx(x:: Array{<:Number}; A=A) = A'*(A*x)

    Lf = real(eigsolve(ATAx, n, 1, :LM, eltype(A))[1][1])
    ϵ  = ϵ₀*Lf
    L  = 1.01*Lf
    λ  = 0.1*norm(A'b, Inf)^2/(2*L) # Pode ser c*norm(A'b, Inf)^2/(2*L), onde 0<c<1
    x₀ = zeros(Float64, n)

    println("Lf, λ: ", Lf, ", ", λ)

    f(x:: Array{<:Number}; A=A, b=b)    = norm(A*x.-b)^2/2
    h(x:: Array{<:Number}; λ=λ)         = λ*norm(x, 0)
    h1(x:: Array{<:Number}; λ1=λ1)      = λ1*norm(x, 1)
    F(x:: Array{<:Number}; f=f, h=h)    = f(x)+h(x)
    F1(x:: Array{<:Number}; f=f, h1=h1) = f(x)+h1(x)
    ∇f(x:: Array{<:Number}; A=A, b=b)   = A'*(A*x.-b)

    #Lₖl1(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n, λ1=λ1) = L, proxhL1(L, x.-∂fx./L, λ1)

    #x_FISTA, histF, histnψ = FISTA(F1, ∇f, Lₖl1, x₀, L, k_max; ϵ=ϵ)
    #println("l1 FISTA: ", 1-norm(x_FISTA, 0)/n, " & ", f(x_FISTA))

    Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = L, proxhL(L, x.-∂fx./L, λ)
    pαₖ(αₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/αₖ, x.-αₖ.*∂fx, λ)
    proxα(α:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/α, x.-α.*∂fx, λ)
    ℘hλg(λₖ:: Number, k:: Int64, y:: Array{<:Number}, ∂fz:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/λₖ, y.-λₖ.*∂fz, λ)
    Tλ(λₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, λ=λ) = proxhL(1/λₖ, x.-λₖ.*∂fx, λ)

    converg_T[i, :], F_hist[i, :] = experiment(f, h, F, ∇f, Lₖ, pαₖ, proxα, ℘hλg, Tλ, x₀, n, L, k_max, ϵ)
end

plt = performance_profile(PlotsBackend(), converg_T, methods, title="Performance profile of convergence")
plot!(plt, dpi=600, legend=:bottomright)
#savefig(plt, "Experiments/Plots/Performance/performance_1_L.png")
#FileIO.save("converg_T_1_L.jld2", "converg_T_1_L", converg_T)

pltF = performance_profile(PlotsBackend(), F_hist, methods, title="Performance profile of best function value")
plot!(pltF, dpi=600, legend=:bottomright)
#savefig(pltF, "Experiments/Plots/Performance/performance_F_1_L.png")
#FileIO.save("F_hist_1_L.jld2", "F_hist_1_L", F_hist)

c_f=FileIO.load("converg_T_full_L.jld2", "converg_T_full_L")
F_f=FileIO.load("F_hist_full_L.jld2", "F_hist_full_L")