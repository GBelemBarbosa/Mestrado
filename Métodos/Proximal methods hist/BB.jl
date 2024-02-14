using Plots
using LaTeXStrings

function BB(Φ:: Function, ∂f:: Function, pαₖ:: Function, x₀:: Array{<:Number}, αmin:: Number, αmax:: Number, M:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    x_=x=x₀
    ∂fx=∂f(x)
    αₖ=1.0
    hist=[Φ(x)]
    
    k=1
    while true
        x_, x=x, pαₖ(αₖ, x, ∂fx)
        push!(hist, Φ(x))

        if norm(∂fx, p)<ϵ || k==k_max
            break
        end
        k+=1

        sₖ=x.-x_
        ∂fx_, ∂fx=∂fx, ∂f(x)
        αₖ=min(αmax, max(αmin, sₖ'*(∂fx.-∂fx_)/sₖ'sₖ))
    end 

    return x, hist
end