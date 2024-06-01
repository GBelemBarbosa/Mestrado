function NSPGHZd(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, η:: Number, γ:: Number, αmin:: Number, αmax:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    start=time()
    x_best=x_=x=x₀
    nsₖ=αₖ=α₀
    sₖ=zeros(Float64, length(x))
    F_best=Fx=c=F(x)
    ∂fx_=∂fx=∂f(x)
    q=1.0
    histnψ=Tuple{Float64, Float64}[]
    pr=0
    histF=[(time()-start, c)]
    
    k=1
    while true
        while true
            x=proxα(αₖ, x_, ∂fx)

            pr+=1

            Fx=F(x)
            sₖ=x.-x_
            nsₖ=sₖ'sₖ

            if Fx+αₖ*γ*nsₖ/2<=c
                break
            end

            αₖ*=ρ

            if αₖ<αmin || isnan(αₖ)
                return x_best, histF, histnψ, pr
            end
        end
        ∂fx_, ∂fx=∂fx, ∂f(x)

        if F_best>Fx
            x_best=x
            F_best=Fx
        end

        t1=time()
        elapsed=t1-start
        nψ=norm(∂fx.-∂fx_.+(x_.-x)./αₖ, p)
        push!(histF, (elapsed, Fx))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1
        
        q_, q=q, η*q+1
        c=(η*q_*c+Fx)/q
        x_=x
        yₖ=∂fx.-∂fx_
        αₖ=nsₖ/(sₖ'yₖ)
        if αₖ>αmax || αₖ<αmin
            αₖ=sqrt(nsₖ/(yₖ'yₖ))
        end
    end 

    return x_best, histF, histnψ, pr
end