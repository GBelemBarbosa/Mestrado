function BB(Φ:: Function, ∂f:: Function, proxhL:: Function, x₀:: Array{<:Number}, α₀:: Number, αmin:: Number, αmax:: Number, k_max:: Int64; ϵ=eps(), p=Inf)
    start=time()
    x_best=x_=x=x₀
    ∂fx_=∂fx=∂f(x)
    Φ_best=Φx=Φ(x₀)
    αₖ=α₀
    histnψ=Tuple{Float64, Float64}[]
    histF=[(time()-start, Φx)]
    
    k=1
    while true
        x_, x=x, proxhL(1/αₖ, x.-αₖ.*∂fx)
        Φx=Φ(x)
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

        sₖ=x.-x_
        ∂fx_, ∂fx=∂fx, ∂f(x)
        αₖ=min(αmax, max(αmin, sₖ'sₖ/(sₖ'*(∂fx.-∂fx_))))
    end 

    return x_best, histF, histnψ
end