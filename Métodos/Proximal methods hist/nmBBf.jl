function nmBBf(Φ:: Function, ∂f:: Function, pαₖ:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, γ:: Number, αmin:: Number, garbage:: Number, M:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    start=time()
    x_best=x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    Φ_best=Φx=Φ(x₀)
    nsₖ=αₖ=α₀
    sₖ=zeros(Float64, length(x))
    last_M=[Φx for i=1:M]
    histnψ=Tuple{Float64, Float64}[]
    pr=0
    histF=[(time()-start, last_M[begin])]
    
    k=1
    while true
        max_M=maximum(last_M)

        while true
            x=pαₖ(αₖ, x_, ∂fx)

            pr+=1

            Φx=Φ(x)
            sₖ=x.-x_
            nsₖ=sₖ'sₖ

            if Φx+αₖ*γ*nsₖ/2<=max_M 
                break
            end

            αₖ*=ρ

            if αₖ<αmin || isnan(αₖ)
                return x_best, histF, histnψ, pr
            end
        end
        ∂fx_, ∂fx=∂fx, ∂f(x)

        if Φ_best>Φx
            x_best=x
            Φ_best=Φx
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

        popfirst!(last_M)
        push!(last_M, Φx)
        x_=x
        αₖ=nsₖ/(sₖ'*(∂fx.-∂fx_))
    end 

    return x_best, histF, histnψ, pr
end