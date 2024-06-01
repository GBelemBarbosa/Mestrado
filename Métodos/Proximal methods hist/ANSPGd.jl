function ANSPGd(F:: Function, ∂f:: Function, prox:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, δ:: Number, αmin:: Number, αmax:: Number, m:: Int64, γ₀:: Number, τ:: Number, β:: Number, γmin:: Number, γmax:: Number, n:: Int64, k_max:: Int64; ϵ=eps(), p=Inf)
    t1=time()
    ∂fx=∂f(x₀)
    start=time()
    maybe=t1-start
    x_best=v=y=z=x_2=x_=x=x₀
    step=γ=γ₀
    F_best=Fv=Fz=Fx=c=F(x)
    ∂fx_=∂fx=∂f(x)
    ∂fy_=∂fy=∂f(y)
    last_n=[Fx for i=1:n]
    last_m=[Fx for i=1:m]
    nzy=t=1.0
    histnψ=Tuple{Float64, Float64}[]
    pr=0
    gr=1
    histF=[(time()-start, c)]
    
    k=1
    while true
        max_n=maximum(last_n)

        while true
            z=prox(γ, y, ∂fy)

            pr+=1

            Fz=F(z)
            nzy=norm(z.-y)^2

            if Fz+β*γ*nzy<=max_n 
                break
            end

            γ*=τ

            if isnan(γ) || γ<γmin
                return x_best, histF, histnψ, pr, gr
            end
        end

        max_m=maximum(last_m)

        if Fz+β*γ*nzy<=max_m
            x=z
            Fx=Fz
            ∂fx=∂fy
            x_2=y
            ∂fx_=∂fy
            step=γ
        else
            if k>1
                xsₖ=x.-y
                nxsₖ=xsₖ'xsₖ
                xyₖ=∂fx.-∂fy_
                α=nxsₖ/(xsₖ'xyₖ)
                if α>αmax || α<αmin
                    α=sqrt(nxsₖ/(xyₖ'xyₖ))
                end
            else
                α=α₀
            end

            start+=maybe
            gr+=1

            while true
                v=prox(α, x, ∂fx)

                pr+=1

                Fv=F(v)

                if Fv+δ*α*norm(v.-x)^2<=max_m
                    break
                end
            
                α*=ρ

                if isnan(α) || α<αmin
                    return x_best, histF, histnψ, pr, gr
                end
            end
            
            if Fv<Fz
                x_2=x
                x=v
                Fx=Fv
                ∂fx_=∂fx
                step=α
            else
                x=z
                Fx=Fz
                x_2=y
                ∂fx_=∂fy
                step=γ
            end
        end
            
        t1=time()
        ∂fx=∂f(x)
        maybe=t1-time()

        if F_best>Fx
            x_best=x
            F_best=Fx
        end

        elapsed=t1-start
        nψ=norm(∂fx.-∂fx_.+(x_2.-x)./step, p)
        push!(histF, (elapsed, Fx))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1
        gr+=1
        
        t_, t=t, (1+sqrt(1+4*t^2))/2
        y_, y=y, x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
        ysₖ=y.-y_
        nysₖ=ysₖ'ysₖ
        ∂fy_, ∂fy=∂fy, ∂f(y)
        popfirst!(last_n)
        push!(last_n, F(y))
        popfirst!(last_m)
        push!(last_m, Fx)
        yyₖ=∂fy.-∂fy_
        γ=nysₖ/(ysₖ'yyₖ)
        if γ>γmax || γ<γmin
            γ=sqrt(nysₖ/(yyₖ'yyₖ))
        end
    end 

    return x_best, histF, histnψ, pr, gr
end