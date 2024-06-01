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

methods = ["NSPG", "-", "NSPGHZ", "η = 0.1", "η = 0.2", "η = 0.5", "η = 0.8", "η = 0.9"]
phantom = count(x->occursin("NSPG", x), methods)

const k_max   = 5000
const ϵ₀      = 10^-5

T_hist  = Array{Float64}(undef, 25, length(methods)-phantom)
F_hist  = Array{Float64}(undef, 25, length(methods)-phantom)
pr_hist = Array{Float64}(undef, 25, length(methods)-phantom)
gr_hist = Array{Float64}(undef, 25, length(methods)-phantom)

# Localização dos dados
lassodir = "../../Data-Lasso/"

for i=1:25
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

    T_hist[i, :], F_hist[i, :], pr_hist[i, :], gr_hist[i, :] = experiment(f, h, F, ∇f, Lₖ, pαₖ, proxα, ℘hλg, Tλ, x₀, n, L, k_max, ϵ)
end
At = findall(x->occursin("NSPG", x), methods)

for at∈At
    for i=at+1:length(methods)
        if occursin("-", methods[i])
            methods[i]=methods[at]*lstrip(methods[i], '-')
        elseif occursin("=", methods[i])
            continue
        else
            break
        end
    end
end
deleteat!(methods, At)

#T_hist  = FileIO.load("Experiments/Plots/Performance/data/T_hist_nmSPG_.jld2", "T_hist_nmSPG")
#F_hist  = FileIO.load("Experiments/Plots/Performance/data/F_hist_nmSPG_.jld2", "F_hist_nmSPG")
#pr_hist = FileIO.load("Experiments/Plots/Performance/data/pr_hist_nmSPG_.jld2", "pr_hist_nmSPG")
#gr_hist = FileIO.load("Experiments/Plots/Performance/data/gr_hist_nmSPG_.jld2", "gr_hist_nmSPG")

pltpr = performance_profile(PlotsBackend(), pr_hist, methods, title="Performance profile of prox calculations")
plot!(pltpr, dpi=600, legend=:bottomright)
savefig(pltpr, "Experiments/Plots/Performance/performance_pr_HZ_eta.png")
FileIO.save("Experiments/Plots/Performance/data/pr_hist_HZ_eta.jld2", "pr_hist_nmSPG", pr_hist)

pltgr = performance_profile(PlotsBackend(), gr_hist, methods, title="Performance profile of gradient evaluations")
plot!(pltgr, dpi=600, legend=:bottomright)
savefig(pltgr, "Experiments/Plots/Performance/performance_gr_HZ_eta.png")
FileIO.save("Experiments/Plots/Performance/data/gr_hist_HZ_eta.jld2", "gr_hist_nmSPG", gr_hist)

pltT = performance_profile(PlotsBackend(), T_hist, methods, title="Performance profile of convergence time")
plot!(pltT, dpi=600, legend=:bottomright)
savefig(pltT, "Experiments/Plots/Performance/performance_HZ_eta.png")
FileIO.save("Experiments/Plots/Performance/data/T_hist_HZ_eta.jld2", "T_hist_nmSPG", T_hist)

pltF = performance_profile(PlotsBackend(), F_hist, methods, title="Performance profile of best function value")
plot!(pltF, dpi=600, legend=:bottomright)
savefig(pltF, "Experiments/Plots/Performance/performance_F_HZ_eta.png")
FileIO.save("Experiments/Plots/Performance/data/F_hist_HZ_eta.jld2", "F_hist_nmSPG", F_hist)
