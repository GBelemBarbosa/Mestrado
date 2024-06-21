function nmBBadp(Φ:: Function, ∂f:: Function, pαₖ:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, γ:: Number, αmin:: Number, αmax:: Number, L:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    start=time()
    x_best=x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    Φ_c=Φ_best=Φx=Φ(x₀)
    Φ_r=Inf
    nsₖ=αₖ=α₀
    sₖ=zeros(Float64, length(x))
    l=pr=0
    histnψ=Tuple{Float64, Float64}[]
    histF=[(time()-start, Φx)]
    
    k=1
    while true
        while true
            x=pαₖ(αₖ, x_, ∂fx)

            pr+=1

            Φx=Φ(x)
            sₖ=x.-x_
            nsₖ=sₖ'sₖ

            if Φx+αₖ*γ*nsₖ/2<=Φ_r 
                break
            end

            αₖ*=ρ

            if isnan(αₖ) || αₖ<αmin 
                return x_best, histF, histnψ, pr
            end
        end
        ∂fx_, ∂fx=∂fx, ∂f(x)

        if Φ_best>Φx
            x_best=x
            Φ_best=Φ_c=Φx
            l=0
        else
            Φ_c=max(Φ_c, Φx)
            l+=1
            
            if l==L
                Φ_r=Φ_c 
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

        x_=x
        αₖ=min(αmax, max(αmin, nsₖ/(sₖ'*(∂fx.-∂fx_))))
    end 

    return x_best, histF, histnψ, pr
end