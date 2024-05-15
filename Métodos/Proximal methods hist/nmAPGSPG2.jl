function nmAPGSPG2(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, α₀:: Number, η:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    t1=time()
    ∂fx_=∂fx=∂f(x₀)
    maybe=t1-time()
    start=time()
    x_best=y=z=x_=x=x₀
    αy=α₀
    F_best=c=F(x)
    q=t=1.0
    histnψ=Tuple{Float64, Float64}[]
    ls=0
    histF=[(time()-start, c)]
    
    k=1
    while true
        ∂fy=∂f(y)
        z=proxα(αy, y, ∂fy)
        
        ls+=1
        
        Fz=F(z)
        ysₖ=z.-y
        nysₖ=ysₖ'ysₖ
        if Fz+δ*nysₖ<=c
            x=z
            Fx=Fz
            ∂fx=∂fy
            x_2=y
            ∂fx_2=∂fy
            αx_=αy
        else
            if k>1
                xsₖ=x.-x_
                αx=xsₖ'xsₖ/(xsₖ'*(∂fx.-∂fx_))
            else
                αx=α₀
            end
            
            v=proxα(αx, x, ∂fx)

            ls+=1
            start+=maybe

            Fv=F(v)
            if Fz>Fv
                x=v
                Fx=Fv
                x_2=x
                ∂fx_2=∂fx
                αx_=αx
            else
                x=z
                Fx=Fz
                x_2=y
                ∂fx_2=∂fy
                αx_=αy
            end
        end
        t1=time()
        ∂fx_, ∂fx=∂fx, ∂f(x)
        maybe=t1-time()

        if F_best>Fx
            x_best=x
            F_best=Fx
        end

        elapsed=t1-start
        nψ=norm(∂fx.-∂fx_2.+(x_2.-x).*αx_, p)
        push!(histF, (elapsed, Fx))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1
        
        t_, t=t, (1+sqrt(1+4*t^2))/2
        q_, q=q, η*q+1
        c=(η*q_*c+Fx)/q
        y=x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
        αy=nysₖ/(ysₖ'*(∂f(z).-∂fy))
    end 

    return x_best, histF, histnψ, ls
end