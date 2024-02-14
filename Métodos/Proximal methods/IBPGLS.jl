using Plots
using LaTeXStrings

function IBPGLS(F:: Function, ∂f:: Function, ℘hλg:: Function, α°ₖ:: Function, β°ₖ:: Function, λ°ₖ:: Function, pₖ:: Function, x₀:: Array{<:Number}, η₁:: Number, η₂:: Number, τ::number, d:: Number, δ:: Number, λₘᵢₙ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    z=x_=x=x₀
    Eδ=Ẽ=F(x)
    αₖ=βₖ=λₖ=0
    
    k=1
    while true
        s=x.-x_
        x_=copy(x)
        αₖ=α°ₖ(k, αₖ)
        βₖ=β°ₖ(k, βₖ)
        λₖ=λ°ₖ(k, λₖ)

        while true
            z=x.+βₖ.*s
            ∂fz=∂f(z)
            x=℘hλg(αₖ, λₖ, k, x_, x, ∂fz) 

            nx=norm(x.-x_)^2
            Eδ=F(x)+δ*nx/(4*λₖ)
            if Eδ<=Ẽ-d*nx/2
                break
            end

            αₖ*=η₁
            βₖ*=η₂ 
            λₖ=max(τ*λₖ, λₘᵢₙ)
        end

        if norm(∂fz, p)<ϵ || k==k_max
            break
        end
        k+=1

        pk=pₖ(k)
        Ẽ=pk*Eδ+(1-pk)*Ẽ
    end 

    return x
end