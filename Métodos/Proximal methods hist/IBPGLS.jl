using Plots
using LaTeXStrings

function IBPGLS(F:: Function, ∂f:: Function, ℘hλg:: Function, α°ₖ:: Function, β°ₖ:: Function, λ°ₖ:: Function, pₖ:: Function, x₀:: Array{<:Number}, η₁:: Number, η₂:: Number, τ:: Number, d:: Number, δ:: Number, λₘᵢₙ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=∂fz=x_=x=x₀
    Eδ=Ẽ=F(x)
    nx=αₖ=βₖ=λₖ=0.0
    hist=[Eδ]
    histnψ=[]
    ls=0
    
    k=1
    while true
        s=x.-x_
        αₖ=α°ₖ(k, αₖ)
        βₖ=β°ₖ(k, βₖ)
        λₖ=λ°ₖ(k, λₖ)
        xₐ=x

        while true
            y=xₐ.+αₖ.*s
            z=xₐ.+βₖ.*s
            ∂fz=∂f(z)
            x=℘hλg(αₖ, λₖ, k, y, ∂fz) 

            ls+=1

            nx=norm(x.-xₐ)^2
            Eδ=F(x)+δ*nx/(4*λₖ)
            if Eδ+d*nx/2<=Ẽ
                break
            end

            αₖ*=η₁
            βₖ*=η₂ 
            λₖ=max(τ*λₖ, λₘᵢₙ)

            if isnan(αₖ) || isnan(βₖ) || min(αₖ, βₖ)<=ϵ
                return x_, hist, histnψ, ls
            end
        end

        push!(hist, Eδ-δ*nx/(4*λₖ))
        nψ=norm(∂f(x).-∂fz.+(y.-x)./λₖ, p)
        push!(nψ, histnψ)

        if nψ<ϵ || k==k_max
            break
        end
        k+=1

        x_=xₐ
        pk=pₖ(k)
        Ẽ=pk*Eδ+(1-pk)*Ẽ
    end 

    return x, hist, histnψ, ls
end