include("../Métodos/Proximal methods hist/BB.jl")
include("../Métodos/Proximal methods hist/nmBB.jl")
include("../Métodos/Proximal methods hist/FISTA.jl")
include("../Métodos/Proximal methods hist/proximal_subgradient.jl")
include("../Métodos/Proximal methods hist/mAPG.jl")
include("../Métodos/Proximal methods hist/nmAPG.jl")
include("../Métodos/Proximal methods hist/nmAPGSPG.jl")
include("../Métodos/Proximal methods hist/nmAPGSPG2.jl")
include("../Métodos/Proximal methods hist/mAPGLS.jl")
include("../Métodos/Proximal methods hist/nmAPGLS.jl")
include("../Métodos/Proximal methods hist/nmAPGLS2.jl")
include("../Métodos/Proximal methods hist/IBPGLS.jl")
include("../Métodos/Proximal methods hist/newAPG_vs.jl")
include("../Métodos/Proximal methods hist/newAPG_vs_x_2.jl")
include("../Métodos/Proximal methods hist/newAPG_vs_y_2.jl")

function experiment(f:: Function, h:: Function, F:: Function, ∇f:: Function, Lₖ:: Function, pαₖ:: Function, proxα:: Function, ℘hλg:: Function, Tλ:: Function, x₀:: Array{<:Number}, n:: Int64, L:: Number, k_max:: Int64, ϵ:: Number)
    converg_T_i = Float64[]
    F_i         = Float64[]

    α₀ = norm(∇f(x₀).-∇f(x₀.+10^-5))/(sqrt(n)*10^-5)
    if "SPG"∈methods
        x_SPG, histF, histnψ = BB(F, ∇f, pαₖ, x₀, α₀, 10^-30, 10^30, k_max; ϵ=ϵ)   
        println("x_SPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_SPG, 0)/n, ", ", f(x_SPG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_SPG))
    end

    # Não sei a razão, mas o primeiro nmBB rodado tem sempre um atraso, então rodo com 1 iteração e descarto
    nmBB(F, ∇f, pαₖ, x₀, α₀, 4, 0.01, 10^-30, 10^30, 1, 1; ϵ=ϵ)

    if "SPGmLS"∈methods
        x_mSPG, histF, histnψ, ls_mSPG = nmBB(F, ∇f, pαₖ, x₀, α₀, 4, 0.01, 10^-30, 10^30, 1, k_max; ϵ=ϵ)
        println("x_SPGmLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_mSPG, 0)/n, ", ", f(x_mSPG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_mSPG))
    end

    if "SPGnmLS"∈methods
        x_nmSPG, histF, histnψ, ls_nmSPG = nmBB(F, ∇f, pαₖ, x₀, α₀, 4, 0.01, 10^-30, 10^30, 5, k_max; ϵ=ϵ)
        println("x_SPGnmLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmSPG, 0)/n, ", ", f(x_nmSPG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmSPG))
    end

    if "mAPG"∈methods
        x_mAPG, histF, histnψ = mAPG(F, ∇f, proxα, x₀, 1/L, 1/L, k_max; ϵ=ϵ) # 2 prox steps per iteration
        println("x_mAPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_mAPG, 0)/n, ", ", f(x_mAPG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_mAPG))
    end

    if "nmAPG"∈methods
        x_nmAPG, histF, histnψ, ls_nmAPG = nmAPG(F, ∇f, proxα, x₀, 1/L, 1/L, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPG, 0)/n, ", ", f(x_nmAPG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPG))
    end

    if "nmAPGSPG"∈methods
        x_nmAPGSPG, histF, histnψ, ls_nmAPGSPG = nmAPGSPG(F, ∇f, proxα, x₀, α₀, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGSPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGSPG, 0)/n, ", ", f(x_nmAPGSPG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGSPG))
    end

    if "nmAPGSPG2"∈methods
        x_nmAPGSPG2, histF, histnψ, ls_nmAPGSPG2 = nmAPGSPG2(F, ∇f, proxα, x₀, α₀, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGSPG2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGSPG2, 0)/n, ", ", f(x_nmAPGSPG2))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGSPG2))
    end

    if "mAPGLS"∈methods
        x_mAPGLS, histF, histnψ, ls_mAPGLS = mAPGLS(F, ∇f, proxα, x₀, α₀, 2/5, 10^-4, k_max; ϵ=ϵ)
        println("x_mAPGLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_mAPGLS, 0)/n, ", ", f(x_mAPGLS))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_mAPGLS))
    end

    if "nmAPGLS"∈methods
        x_nmAPGLS, histF, histnψ, ls_nmAPGLS = nmAPGLS(F, ∇f, proxα, x₀, α₀, 2/5, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGLS, 0)/n, ", ", f(x_nmAPGLS))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGLS))
    end

    if "nmAPGLS2"∈methods
        x_nmAPGLS2, histF, histnψ, ls_nmAPGLS2 = nmAPGLS2(F, ∇f, proxα, x₀, α₀, 2/5, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGLS2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGLS2, 0)/n, ", ", f(x_nmAPGLS2))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGLS2))
    end

    if "PG"∈methods
        x_PG, histF, histnψ = PG(F, ∇f, Lₖ, x₀, L, k_max; ϵ=ϵ)
        println("x_PG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_PG, 0)/n, ", ", f(x_PG))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_PG))
    end

    if "FISTA"∈methods
        x_FISTA, histF, histnψ = FISTA(F, ∇f, Lₖ, x₀, L, k_max; ϵ=ϵ)
        println("x_FISTA:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_FISTA, 0)/n, ", ", f(x_FISTA))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_FISTA))
    end

    α°ₖ(k:: Int64, αₖ:: Number)=9
    β°ₖ(k:: Int64, βₖ:: Number)=0.9
    λ°ₖ(k:: Int64, λₖ:: Number)=0.5 
    λₘᵢₙ=0.98*(1-0.1)/L
    d=((1-0.1)/λₘᵢₙ-L)/2
    pₖ(k:: Int64)=1
    if "IBPGmLS"∈methods
        x_IBPGmLS, histF, histnψ, ls_IBPGmLS = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGmLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGmLS, 0)/n, ", ", f(x_IBPGmLS))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGmLS))
    end

    λ°ₖL(k:: Int64, λₖ:: Number; L=L)=0.99*(1-0.1)/L
    if "IBPGmLSL"∈methods
        x_IBPGmLSL, histF, histnψ, ls_IBPGmLSL = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖL, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGmLSL:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGmLSL, 0)/n, ", ", f(x_IBPGmLSL))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGmLSL))
    end

    pₖ(k:: Int64)=0.7
    if "IBPGnmLS"∈methods
        x_IBPGnmLS, histF, histnψ, ls_IBPGnmLS = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGnmLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGnmLS, 0)/n, ", ", f(x_IBPGnmLS))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGnmLS))
    end

    if "IBPGnmLSL"∈methods
        x_IBPGnmLSL, histF, histnψ, ls_IBPGnmLSL = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖL, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGnmLSL:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGnmLSL, 0)/n, ", ", f(x_IBPGnmLSL))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGnmLSL))
    end

    function γ(k:: Int64, t_:: Number)
        t=(1+sqrt(1+4*t_^2))/2

        return t, (t_-1)/t
    end
    Q(k:: Int64)=0.99^k
    E(k:: Int64)=k^-1.1
    λ₁=1/α₀
    if "newAPG_vs"∈methods
        x_newAPG_vs, histF, histnψ, ls_newAPG_vs = newAPG_vs(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
        println("x_newAPG_vs:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_newAPG_vs, 0)/n, ", ", f(x_newAPG_vs))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_newAPG_vs))
    end

    if "newAPG_vs_x_2"∈methods
        x_newAPG_vs_x_2, histF, histnψ, ls_newAPG_vs_x_2 = newAPG_vs_x_2(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
        println("x_newAPG_vs_x_2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_newAPG_vs_x_2, 0)/n, ", ", f(x_newAPG_vs_x_2))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_newAPG_vs_x_2))
    end

    if "newAPG_vs_y_2"∈methods
        x_newAPG_vs_y_2, histF, histnψ, ls_newAPG_vs_y_2 = newAPG_vs_y_2(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
        println("x_newAPG_vs_y_2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_newAPG_vs_y_2, 0)/n, ", ", f(x_newAPG_vs_y_2))
        push!(converg_T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_newAPG_vs_y_2))
    end

    return converg_T_i, F_i
end