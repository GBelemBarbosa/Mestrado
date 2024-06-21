include("../Métodos/Proximal methods hist/BB.jl")
include("../Métodos/Proximal methods hist/nmBB.jl")
include("../Métodos/Proximal methods hist/nmBBf.jl")
include("../Métodos/Proximal methods hist/nmBBd.jl")
include("../Métodos/Proximal methods hist/nmBBalt.jl")
include("../Métodos/Proximal methods hist/nmBBadp.jl")
include("../Métodos/Proximal methods hist/nmBBadp2.jl")
include("../Métodos/Proximal methods hist/CBB.jl")
include("../Métodos/Proximal methods hist/NPGLSHZd.jl")
include("../Métodos/Proximal methods hist/ANSPG.jl")
include("../Métodos/Proximal methods hist/ANSPGf.jl")
include("../Métodos/Proximal methods hist/ANSPGd.jl")
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
    T_i  = Float64[]
    F_i  = Float64[]
    pr_i = Float64[]
    gr_i = Float64[]

    α₀ = (sqrt(n)*10^-5)/norm(∇f(x₀).-∇f(x₀.+10^-5))
    println("α₀ = ", α₀)
    if "SPG"∈methods
        x_SPG, histF, histnψ = BB(F, ∇f, pαₖ, x₀, α₀, α₀/1000, 10^30, k_max; ϵ=ϵ)   
        println("x_SPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_SPG, 0)/n, ", ", f(x_SPG))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_SPG))
    end

    # Não sei a razão, mas o primeiro nmBB rodado tem sempre um atraso, então rodo com 1 iteração e descarto
    nmBB(F, ∇f, pαₖ, x₀, α₀, 1/2, 0.01, α₀/1000, 10^30, 1, 1; ϵ=ϵ)

    if "NSPG"∈methods
        at = findfirst(x->x=="NSPG", methods)
        
        for i=at+1:length(methods)
            m = 5
            ρ = 1/2
            γ = 0.01

            if occursin("m = ", methods[i]) 
                m = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("ρ = ", methods[i])
                ρ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("γ = ", methods[i])
                γ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif !occursin("-", methods[i])
                break
            end
            if occursin("f", methods[i])
                NSPG = nmBBf
            elseif occursin("d", methods[i])
                NSPG = nmBBd
            elseif occursin("a", methods[i])
                NSPG = nmBBalt
            else
                NSPG = nmBB
            end

            x_NSPG, histF, histnψ, pr_NSPG = NSPG(F, ∇f, pαₖ, x₀, α₀, ρ, γ, α₀/1000, 10^30, m, k_max; ϵ=ϵ)
            println("x_NSPG ("*methods[i]*")")
            println("pr, gr, nψ[end], sparsity, f_best: ", pr_NSPG+Inf*(histnψ[end][2]>=ϵ), ", ", length(histnψ)+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_NSPG, 0)/n, ", ", f(x_NSPG))
            push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
            push!(F_i, F(x_NSPG))
            push!(pr_i, pr_NSPG+Inf*(histnψ[end][2]>=ϵ))
            push!(gr_i, length(histnψ)+Inf*(histnψ[end][2]>=ϵ))
        end
    end

    if "NSPGadp"∈methods
        at = findfirst(x->x=="NSPGadp", methods)
        
        for i=at+1:length(methods)
            L = 5
            ρ = 1/2
            γ = 0.01

            if occursin("L = ", methods[i]) 
                L = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("ρ = ", methods[i])
                ρ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("γ = ", methods[i])
                γ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif !occursin("-", methods[i])
                break
            end

            x_NSPGadp, histF, histnψ, pr_NSPGadp = nmBBadp(F, ∇f, pαₖ, x₀, α₀, ρ, γ, α₀/1000, 10^30, L, k_max; ϵ=ϵ)
            println("x_NSPGadp ("*methods[i]*")")
            println("pr, gr, nψ[end], sparsity, f_best: ", pr_NSPGadp+Inf*(histnψ[end][2]>=ϵ), ", ", length(histnψ)+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_NSPGadp, 0)/n, ", ", f(x_NSPGadp))
            push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
            push!(F_i, F(x_NSPGadp))
            push!(pr_i, pr_NSPGadp+Inf*(histnψ[end][2]>=ϵ))
            push!(gr_i, length(histnψ)+Inf*(histnψ[end][2]>=ϵ))
        end
    end

    if "NSPGadp2"∈methods
        at = findfirst(x->x=="NSPGadp2", methods)
        
        for i=at+1:length(methods)
            m = 6
            L = 5
            P = 4*m
            ρ = 1/2
            γ = 0.01

            if occursin("M = ", methods[i]) 
                m = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("L = ", methods[i]) 
                L = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("P = ", methods[i]) 
                P = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("ρ = ", methods[i])
                ρ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("γ = ", methods[i])
                γ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif !occursin("-", methods[i])
                break
            end
            γ₁ = m/L
            γ₂ = P/m

            x_NSPGadp2, histF, histnψ, pr_NSPGadp2 = nmBBadp2(F, ∇f, pαₖ, x₀, α₀, ρ, γ, α₀/1000, 10^30, m, L, P, γ₁, γ₂, k_max; ϵ=ϵ)
            println("x_NSPGadp2 ("*methods[i]*")")
            println("pr, gr, nψ[end], sparsity, f_best: ", pr_NSPGadp2+Inf*(histnψ[end][2]>=ϵ), ", ", length(histnψ)+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_NSPGadp2, 0)/n, ", ", f(x_NSPGadp2))
            push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
            push!(F_i, F(x_NSPGadp2))
            push!(pr_i, pr_NSPGadp2+Inf*(histnψ[end][2]>=ϵ))
            push!(gr_i, length(histnψ)+Inf*(histnψ[end][2]>=ϵ))
        end
    end

    if "CBB"∈methods
        at = findfirst(x->x=="CBB", methods)
        
        for i=at+1:length(methods)
            m = 5
            L = 5
            P = 5
            N = 4
            ρ = 1/2
            γ = 0.01

            if occursin("m = ", methods[i]) 
                m = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("L = ", methods[i]) 
                L = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("P = ", methods[i]) 
                P = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("N = ", methods[i]) 
                N = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("ρ = ", methods[i])
                ρ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("γ = ", methods[i])
                γ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("θ = ", methods[i])
                θ  = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif !occursin("-", methods[i])
                break
            end
            γ₁ = m/L
            γ₂ = P/m

            x_CBB, histF, histnψ, pr_CBB = CBB(F, ∇f, pαₖ, x₀, α₀, ρ, γ, α₀/1000, 10^30, m, L, P, N, γ₁, γ₂, θ, k_max; ϵ=ϵ)
            println("x_CBB ("*methods[i]*")")
            println("pr, gr, nψ[end], sparsity, f_best: ", pr_CBB+Inf*(histnψ[end][2]>=ϵ), ", ", length(histnψ)+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_CBB, 0)/n, ", ", f(x_CBB))
            push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
            push!(F_i, F(x_CBB))
            push!(pr_i, pr_CBB+Inf*(histnψ[end][2]>=ϵ))
            push!(gr_i, length(histnψ)+Inf*(histnψ[end][2]>=ϵ))
        end
    end

    if "NSPGHZ"∈methods
        at = findfirst(x->x=="NSPGHZ", methods)
        
        for i=at+1:length(methods)
            ρ = 1/2
            γ = 0.01
            η = 0.5

            if occursin("ρ = ", methods[i])
                ρ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("γ = ", methods[i])
                γ = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif occursin("η = ", methods[i])
                η = parse(Float64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif !occursin("-", methods[i])
                break
            end
            if occursin("f", methods[i])
                NSPG = NSPGHZf
            elseif occursin("a", methods[i])
                NSPG = NSPGHZa
            else
                NSPG = NSPGHZd
            end

            x_NSPGHZ, histF, histnψ, pr_NSPGHZ = NSPG(F, ∇f, pαₖ, x₀, α₀, ρ, η, γ, α₀/1000, 10^30, k_max; ϵ=ϵ)
            println("x_NSPGHZ ("*methods[i]*")")
            println("pr, gr, nψ[end], sparsity, f_best: ", pr_NSPGHZ+Inf*(histnψ[end][2]>=ϵ), ", ", length(histnψ)+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_NSPGHZ, 0)/n, ", ", f(x_NSPGHZ))
            push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
            push!(F_i, F(x_NSPGHZ))
            push!(pr_i, pr_NSPGHZ+Inf*(histnψ[end][2]>=ϵ))
            push!(gr_i, length(histnψ)+Inf*(histnψ[end][2]>=ϵ))
        end
    end

    if "ANSPG"∈methods
        at = findfirst(x->x=="ANSPG", methods)
        
        for i=at+1:length(methods)
            m = n₂ = 5
            ρ = τ  = 1/2
            δ = β  = 0.01

            if occursin("n =", methods[i])
                n₂     = parse(Int64, methods[i][findfirst("=", methods[i])[1]+2:end])
            elseif !occursin("-", methods[i])
                break
            end
            if occursin("f", methods[i])
                ANSPGa = ANSPGf
            elseif occursin("d", methods[i])
                ANSPGa = ANSPGd
            else
                ANSPGa = ANSPG
            end

            x_ANSPG, histF, histnψ, pr_ANSPG, gr_ANSPG = ANSPGa(F, ∇f, pαₖ, x₀, α₀, ρ, δ, α₀/1000, 10^30, m, α₀, τ, β, α₀/1000, 10^30, n₂, k_max; ϵ=ϵ)
            println("x_ANSPG ("*methods[i]*")")
            println("pr, gr, nψ[end], sparsity, f_best: ", pr_ANSPG+Inf*(histnψ[end][2]>=ϵ), ", ", gr_ANSPG+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_ANSPG, 0)/n, ", ", f(x_ANSPG))
            push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
            push!(F_i, F(x_ANSPG))
            push!(pr_i, pr_ANSPG+Inf*(histnψ[end][2]>=ϵ))
            push!(gr_i, gr_ANSPG+Inf*(histnψ[end][2]>=ϵ))
        end
    end

    if "nmAPGLS"∈methods
        x_nmAPGLS, histF, histnψ, pr_nmAPGLS, gr_nmAPGLS = nmAPGLS(F, ∇f, proxα, x₀, α₀, 2/5, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGLS:")
        println("pr, gr, nψ[end], sparsity, f_best: ", pr_nmAPGLS+Inf*(histnψ[end][2]>=ϵ), ", ", gr_nmAPGLS+Inf*(histnψ[end][2]>=ϵ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGLS, 0)/n, ", ", f(x_nmAPGLS))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGLS))
        push!(pr_i, pr_nmAPGLS+Inf*(histnψ[end][2]>=ϵ))
        push!(gr_i, gr_nmAPGLS+Inf*(histnψ[end][2]>=ϵ))
    end

    if "mAPG"∈methods
        x_mAPG, histF, histnψ = mAPG(F, ∇f, proxα, x₀, 1/L, 1/L, k_max; ϵ=ϵ) # 2 prox steps per iteration
        println("x_mAPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_mAPG, 0)/n, ", ", f(x_mAPG))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_mAPG))
    end

    if "nmAPG"∈methods
        x_nmAPG, histF, histnψ, pr_nmAPG = nmAPG(F, ∇f, proxα, x₀, 1/L, 1/L, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPG, 0)/n, ", ", f(x_nmAPG))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPG))
    end

    if "nmAPGSPG"∈methods
        x_nmAPGSPG, histF, histnψ, pr_nmAPGSPG = nmAPGSPG(F, ∇f, proxα, x₀, α₀, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGSPG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGSPG, 0)/n, ", ", f(x_nmAPGSPG))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGSPG))
    end

    if "nmAPGSPG2"∈methods
        x_nmAPGSPG2, histF, histnψ, pr_nmAPGSPG2 = nmAPGSPG2(F, ∇f, proxα, x₀, α₀, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGSPG2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGSPG2, 0)/n, ", ", f(x_nmAPGSPG2))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGSPG2))
    end

    if "mAPGLS"∈methods
        x_mAPGLS, histF, histnψ, pr_mAPGLS = mAPGLS(F, ∇f, proxα, x₀, α₀, 2/5, 10^-4, k_max; ϵ=ϵ)
        println("x_mAPGLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_mAPGLS, 0)/n, ", ", f(x_mAPGLS))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_mAPGLS))
    end

    if "nmAPGLS2"∈methods
        x_nmAPGLS2, histF, histnψ, pr_nmAPGLS2 = nmAPGLS2(F, ∇f, proxα, x₀, α₀, 2/5, 0.8, 10^-4, k_max; ϵ=ϵ)
        println("x_nmAPGLS2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_nmAPGLS2, 0)/n, ", ", f(x_nmAPGLS2))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_nmAPGLS2))
    end

    if "PG"∈methods
        x_PG, histF, histnψ = PG(F, ∇f, Lₖ, x₀, L, k_max; ϵ=ϵ)
        println("x_PG:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_PG, 0)/n, ", ", f(x_PG))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_PG))
    end

    if "FISTA"∈methods
        x_FISTA, histF, histnψ = FISTA(F, ∇f, Lₖ, x₀, L, k_max; ϵ=ϵ)
        println("x_FISTA:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_FISTA, 0)/n, ", ", f(x_FISTA))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_FISTA))
    end

    α°ₖ(k:: Int64, αₖ:: Number)=9
    β°ₖ(k:: Int64, βₖ:: Number)=0.9
    λ°ₖ(k:: Int64, λₖ:: Number)=0.5 
    λₘᵢₙ=0.98*(1-0.1)/L
    d=((1-0.1)/λₘᵢₙ-L)/2
    pₖ(k:: Int64)=1
    if "IBPGmLS"∈methods
        x_IBPGmLS, histF, histnψ, pr_IBPGmLS = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGmLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGmLS, 0)/n, ", ", f(x_IBPGmLS))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGmLS))
    end

    λ°ₖL(k:: Int64, λₖ:: Number; L=L)=0.99*(1-0.1)/L
    if "IBPGmLSL"∈methods
        x_IBPGmLSL, histF, histnψ, pr_IBPGmLSL = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖL, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGmLSL:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGmLSL, 0)/n, ", ", f(x_IBPGmLSL))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGmLSL))
    end

    pₖ(k:: Int64)=0.7
    if "IBPGnmLS"∈methods
        x_IBPGnmLS, histF, histnψ, pr_IBPGnmLS = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGnmLS:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGnmLS, 0)/n, ", ", f(x_IBPGnmLS))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGnmLS))
    end

    if "IBPGnmLSL"∈methods
        x_IBPGnmLSL, histF, histnψ, pr_IBPGnmLSL = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖL, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
        println("x_IBPGnmLSL:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_IBPGnmLSL, 0)/n, ", ", f(x_IBPGnmLSL))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_IBPGnmLSL))
    end

    function γ(k:: Int64, tₖ₋₁:: Number)
        tₖ=(1+sqrt(1+4*tₖ₋₁^2))/2

        return tₖ, (tₖ₋₁-1)/tₖ
    end
    Q(k:: Int64)=0.99^k
    E(k:: Int64)=k^-1.1
    λ₁ = α₀
    if "newAPG_vs"∈methods
        x_newAPG_vs, histF, histnψ, pr_newAPG_vs = newAPG_vs(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
        println("x_newAPG_vs:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_newAPG_vs, 0)/n, ", ", f(x_newAPG_vs))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_newAPG_vs))
    end

    if "newAPG_vs_x_2"∈methods
        x_newAPG_vs_x_2, histF, histnψ, pr_newAPG_vs_x_2 = newAPG_vs_x_2(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
        println("x_newAPG_vs_x_2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_newAPG_vs_x_2, 0)/n, ", ", f(x_newAPG_vs_x_2))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_newAPG_vs_x_2))
    end

    if "newAPG_vs_y_2"∈methods
        x_newAPG_vs_y_2, histF, histnψ, pr_newAPG_vs_y_2 = newAPG_vs_y_2(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
        println("x_newAPG_vs_y_2:")
        println("k_end, nψ[end], sparsity, f_best: ", length(histnψ), ", ", histnψ[end][2], ", ", 1-norm(x_newAPG_vs_y_2, 0)/n, ", ", f(x_newAPG_vs_y_2))
        push!(T_i, histnψ[end][1]+Inf*(histnψ[end][2]>=ϵ))
        push!(F_i, F(x_newAPG_vs_y_2))
    end

    return T_i, F_i, pr_i, gr_i
end