using Plots
using LaTeXStrings

function BB(Φ:: Function, ∂f:: Function, dαₖ:: Function, x₀:: Array{<:Number}, ρ:: Number, γ:: Number, αmin:: Number, αmax:: Number, M:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    x=x₀
    ∂fx=∂f(x)
    αₖ=1.0
    sₖ=zeros(Float64, length(x))
    last_M=[Φ(x) for i=1:M]
    
    k=1
    while true
        dₖ=dαₖ(αₖ, x, ∂fx) #pαₖ(αₖ, x, ∂fx).-x
        λ=1.0
        max_M=maximum(last_M)
        nsₖ=αₖ*γ*norm(sₖ)^2/2
        x_, x=x, x.+λ.*dₖ
        Φx=Φ(x)

        while Φx+λ*nsₖ>max_M
            λ.*=ρ
            x=x_.+λ.*dₖ
            Φx=Φ(x)
        end

        if norm(∂fx, p)<ϵ || k==k_max
            break
        end
        k+=1

        popat!(last_M, M)
        push!(last_M, Φ(x))
        sₖ=x.-x_
        ∂fx_, ∂fx=∂fx, ∂f(x)
        αₖ=min(αmax, max(αmin, sₖ'(∂fx.-∂fx_)/sₖ'sₖ))
    end 

    return x
end