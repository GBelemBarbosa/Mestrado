using Plots
using LaTeXStrings

function BB(Φ:: Function, ∂f:: Function, pαₖ:: Function, x₀:: Array{<:Number}, αmin:: Number, αmax:: Number, M:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    x_best=x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    Φ_best=Φ(x)
    αₖ=1.0
    
    k=1
    while true
        x_, x=x, pαₖ(αₖ, x, ∂fx)
        Φx=Φ(x)
        if Φ_best>Φx
            x_best=x
            Φ_best=Φx
        end

        if norm(∂fx.-∂fx_.+(x_.-x).*αₖ, p)<ϵ || k==k_max
            break
        end
        k+=1

        sₖ=x.-x_
        ∂fx_, ∂fx=∂fx, ∂f(x)
        αₖ=min(αmax, max(αmin, sₖ'(∂fx.-∂fx_)/sₖ'sₖ))
    end 

    return x_best
end