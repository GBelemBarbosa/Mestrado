function nmBB(Φ:: Function, ∂f:: Function, pαₖ:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, γ:: Number, αmin:: Number, αmax:: Number, M:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    nsₖ=Φx=αₖ=α₀
    sₖ=zeros(Float64, length(x))
    last_M=[Φ(x) for i=1:M]
    
    k=1
    while true
        max_M=maximum(last_M)

        while true
            x=pαₖ(αₖ, x_, ∂fx)
            Φx=Φ(x)
            sₖ=x.-x_
            nsₖ=sₖ'sₖ

            if Φx+αₖ*γ*nsₖ/2<=max_M
                break
            end

            αₖ.*=ρ

            if αₖ>=αmax || isnan(αₖ)
                return x_, hist
            end
        end

        if norm(∂fx.-∂fx_.+(x_.-x).*αₖ, p)<ϵ || k==k_max
            break
        end
        k+=1

        popfirst!(last_M)
        push!(last_M, Φx)
        sₖ=x.-x_
        x_=x
        ∂fx_, ∂fx=∂fx, ∂f(x)
        αₖ=min(αmax, max(αmin, nsₖ/(sₖ'(∂fx.-∂fx_))))
    end 

    return x
end