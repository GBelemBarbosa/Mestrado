using Random
using StatsBase: sample
using LinearAlgebra
using Plots
using LaTeXStrings

include("Métodos/Proximal methods hist/BB.jl")
include("Métodos/Proximal methods hist/nmBB.jl")
include("Métodos/Proximal methods hist/FISTA.jl")
include("Métodos/Proximal methods hist/mAPG.jl")
include("Métodos/Proximal methods hist/nmAPG.jl")
include("Métodos/Proximal methods hist/IBPGLS.jl")
include("Métodos/Proximal methods hist/newAPG_vs.jl")
include("Métodos/Proximal methods hist/newAPG_vs_x_2.jl")
include("Métodos/Proximal methods hist/newAPG_vs_y_2.jl")
include("Group sparsity/group_sparse_functions.jl")

Random.seed!(1)

m     = [250]
n     = 10 .*m
plots = []
k_max = 100

for i=eachindex(m)
    xₒₚₜ = zeros(Float64, n[i]) 
    xₒₚₜ[sample(1:n[i], Int(m[i]/10); replace=false)] .= 1
    global A   = randn(Float64, (m[i], n[i]))
    global b   = A*xₒₚₜ.+0.01.*randn(Float64, m[i])

    Lf = eigmax(A'A)
    L  = 1.01*Lf
    λ  = Lf/4
    print(λ)

    f(x:: Array{<:Number}; A=A, b=b)  = norm(A*x.-b)^2/2
    h(x:: Array{<:Number}; λ=λ)       = λ*norm(x, 0)
    F(x:: Array{<:Number}; f=f, h=h)  = f(x)+h(x)
    ∇f(x:: Array{<:Number}; A=A, b=b) = A'*(A*x.-b)

    Lₖ(L:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n[i], λ=λ) = L, proxhL(L, x.-∂fx./L, λ, n)
    pαₖ(αₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n[i], λ=λ) = proxhL(αₖ, x.-∂fx./αₖ, λ, n)
    proxα(α:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n[i], λ=λ) = proxhL(α, x.-∂fx./α, λ, n)
    ℘hλg(αₖ:: Number, λₖ:: Number, x_:: Array{<:Number}, x:: Array{<:Number}, ∂fz:: Array{<:Number}; proxhL=proxhL, n=n[i], λ=λ) = proxhL(1/λₖ, x.+αₖ.*(x.-x_).-λₖ.*∂fz, λ, n, )
    Tλ(λₖ:: Number, x:: Array{<:Number}, ∂fx:: Array{<:Number}; proxhL=proxhL, n=n[i], λ=λ) = proxhL(1/λₖ, x.-λₖ.*∂fx, λ, n)

    x₀ = randn(n[i])
    p  = plot(title=L"l_0"*" penalty", xlabel="Iteration number", ylabel=L"F(x_k)")

    global x_BB, hist = BB(F, ∇f, pαₖ, x₀, 10^-30, 10^30, 5, k_max)
    plot!(eachindex(hist), hist, marker=:circle, label="BB")

    global x_nmBB, hist = nmBB(F, ∇f, pαₖ, x₀, 4, 0.01, 10^-30, 10^18, 5, k_max)
    plot!(eachindex(hist), hist, marker=:circle, label="nmBB")

    global x_mAPG, hist = mAPG(F, ∇f, proxα, x₀, L, L, k_max)
    plot!(eachindex(hist), hist, marker=:diamond, label="mAPG")

    global x_nmAPG, hist = nmAPG(F, ∇f, proxα, x₀, L, L, 0.8, 10^-4, k_max)
    plot!(eachindex(hist), hist, marker=:diamond, label="nmAPG")

    global x_FISTA, hist = FISTA(F, ∇f, Lₖ, x₀, L, k_max)
    plot!(eachindex(hist), hist, linestyle=:dash, label="FISTA")

    α°ₖ(k:: Int64, αₖ:: Number)=9
    β°ₖ(k:: Int64, βₖ:: Number)=0.9
    λ°ₖ(k:: Int64, λₖ:: Number)=0.5 # Try 0.99*(1-0.1)/Lf
    λₘᵢₙ=0.98*(1-0.1)/Lf
    d=((1-0.1)/λₘᵢₙ-Lf)/2
    pₖ(k:: Int64)=1
    global x_IBPGmLS, hist = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max)
    plot!(eachindex(hist), hist, linestyle=:dash, label="IBPGmLS")

    pₖ(k:: Int64)=0.7
    global x_IBPGnmLS, hist = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max)
    plot!(eachindex(hist), hist, linestyle=:dash, label="IBPGnmLS")

    # marker=:utriangle
    # marker=:x
    # linestyle=:dash

    plot!(ylims=(0, 1.1*hist[begin]))
    
    push!(plots, p)
end