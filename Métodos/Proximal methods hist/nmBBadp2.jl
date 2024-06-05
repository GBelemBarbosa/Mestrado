function nmBBadp2(Φ:: Function, ∂f:: Function, pαₖ:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, γ:: Number, αmin:: Number, αmax:: Number, M:: Int64, L:: Int64, P:: Int64, γ₁:: Number, γ₂:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    start=time()
    x_best=x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    Φ_max=Φ_c=Φ_best=Φx=Φ(x₀)
    Φ_r=Inf
    nsₖ=αₖ=α₀
    sₖ=zeros(Float64, length(x))
    last_M=[Φx for i=1:M]
    histnψ=Tuple{Float64, Float64}[]
    p=l=pr=0
    histF=[(time()-start, last_M[begin])]
    
    k=1
    while true
        x=pαₖ(αₖ, x_, ∂fx)

        pr+=1

        Φx=Φ(x)
        sₖ=x.-x_
        nsₖ=sₖ'sₖ

        if Φx+αₖ*γ*nsₖ/2<=Φ_r
            p+=1
        else
            p=0
            
            while true
                αₖ*=ρ

                if isnan(αₖ) || αₖ<αmin
                    return x_best, histF, histnψ, pr
                end

                x=pαₖ(αₖ, x_, ∂fx)

                pr+=1

                Φx=Φ(x)
                sₖ=x.-x_
                nsₖ=sₖ'sₖ

                if Φx+αₖ*γ*nsₖ/2<=min(Φ_r, Φ_max) 
                    break
                end
            end
        end
        ∂fx_, ∂fx=∂fx, ∂f(x)
        popfirst!(last_M)
        push!(last_M, Φx)
        Φ_max=maximum(last_M)

        if Φ_best>Φx
            x_best=x
            Φ_best=Φ_c=Φx
            l=0
        else
            Φ_c=max(Φ_c, Φx)
            l+=1
            
            if l==L
                if (Φ_max-Φ_best)/(Φ_c-Φ_best)>γ₁
                    Φ_r=Φ_c 
                else
                    Φ_r=Φ_max
                end
                Φ_c=Φx
                l=0
            end
        end

        t1=time()
        elapsed=t1-start
        nψ=norm(∂fx.-∂fx_.+(x_.-x)./αₖ, p)
        push!(histF, (elapsed, Φx))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1

        if p>=P
            if Φ_max!=Φx && (Φ_r-Φx)/(Φ_max-Φx)>=γ₂
                Φ_r=Φ_max 
            end
        end

        x_=x
        αₖ=min(αmax, max(αmin, nsₖ/(sₖ'*(∂fx.-∂fx_))))
    end 

    return x_best, histF, histnψ, pr
end