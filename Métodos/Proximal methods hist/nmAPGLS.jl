function nmAPGLS(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, α₀:: Number, ρ:: Number, η:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    maybe=0
    start=time()
    x_best=v=y=z=x_2=x_=x=x₀
    αx_=αy=α₀
    F_best=Fv=Fz=Fx=c=F(x)
    ∂fy=∂fx_=∂fx=∂f(x)
    q=t=1.0
    histnψ=Tuple{Float64, Float64}[]
    pr=0
    gr=1
    histF=[(time()-start, c)]
    
    k=1
    while true
        Fy=F(y)

        while true
            z=proxα(αy, y, ∂fy)

            pr+=1

            Fz=F(z)
            nzy=norm(z.-y)^2

            if Fz+δ*nzy<=c
                x=z
                Fx=Fz
                ∂fx=∂fy
                x_2=y
                ∂fx_=∂fy
                αx_=αy
                break
            elseif Fz+δ*nzy<=Fy 
                if k>1
                    xsₖ=x.-y
                    αx=xsₖ'xsₖ/(xsₖ'*(∂fx.-∂fy)) #max/min?
                    
                    gr+=1
                else
                    αx=α₀
                end

                start+=maybe

                while true
                    v=proxα(αx, x, ∂fx)

                    pr+=1

                    Fv=F(v)

                    if Fv+δ*norm(v.-x)^2<=c
                        break
                    end

                    αx*=ρ

                    if isnan(αx) || αx<10^-17
                        return x_best, histF, histnψ, pr, gr
                    end
                end
                if Fz>Fv
                    x_2=x
                    x=v
                    Fx=Fv
                    ∂fx_=∂fx
                    αx_=αx
                else
                    x=z
                    Fx=Fz
                    x_2=y
                    ∂fx_=∂fy
                    αx_=αy
                end
                break
            end

            αy*=ρ

            if isnan(αy) || αy<10^-17
                return x_best, histF, histnψ, pr, gr
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
        nψ=norm(∂fx.-∂fx_.+(x_2.-x)./αx_, p)
        push!(histF, (elapsed, Fx))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1
        gr+=1
        
        t_, t=t, (1+sqrt(1+4*t^2))/2
        q_, q=q, η*q+1
        c=(η*q_*c+Fx)/q
        y_, y=y, x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
        ysₖ=y.-y_
        ∂fy_, ∂fy=∂fy, ∂f(y)
        αy=ysₖ'ysₖ/(ysₖ'*(∂fy.-∂fy_))
    end 

    return x_best, histF, histnψ, pr, gr
end