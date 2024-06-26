function nmBB(Φ:: Function, ∂f:: Function, proxhL:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, γ:: Number, αmin:: Number, αmax:: Number, M:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    start=time()
    x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    nsₖ=Φx=αₖ=α₀
    sₖ=zeros(Float64, length(x))
    last_M=[Φ(x) for i=1:M]
    histnψ=Tuple{Float64, Float64}[]
    ls=0
    histF=[(time()-start, last_M[begin])]
    
    k=1
    while true
        max_M=maximum(last_M)

        while true
            x=proxhL(1/αₖ, x_.-αₖ.*∂fx)

            ls+=1

            Φx=Φ(x)
            sₖ=x.-x_
            nsₖ=sₖ'sₖ

            if Φx+αₖ*γ*nsₖ/2<=max_M 
                break
            end

            αₖ*=ρ

            if αₖ>=αmax || isnan(αₖ)
                return x_, histF, histnψ, ls
            end
        end
        ∂fx_, ∂fx=∂fx, ∂f(x)

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
        αₖ=min(αmax, max(αmin, nsₖ/(sₖ'*(∂fx.-∂fx_))))
    end 

    return x, histF, histnψ, ls
end