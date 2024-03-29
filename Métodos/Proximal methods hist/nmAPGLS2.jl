using Plots
using LaTeXStrings

function nmAPGLS2(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, ρ:: Number, η:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    t1=time()
    ∂fx_=∂fx=∂f(x₀)
    start=time()
    maybe=t1-start
    v=ysₖ=y=z=x_2=x_=x=x₀
    Fv=Fz=Fx=c=F(x)
    q=nysₖ=αx_=αy=t=1.0
    histnψ=Tuple{Float64, Float64}[]
    ls=0
    histF=[(time()-start, c)]
    
    k=1
    while true
        Fy=F(y)
        ∂fy=∂f(y)

        while true
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
                ∂fx_=∂fy
                αx_=αy
                break
            elseif Fz+δ*nysₖ<=Fy 
                if k>1
                    xsₖ=x.-y
                    αx=xsₖ'xsₖ/(xsₖ'*(∂fx.-∂fy))
                else
                    αx=1.0
                end

                while true
                    v=proxα(αx, x, ∂fx)

                    ls+=1
                    start+=maybe
                    
                    Fv=F(v)

                    if Fv+δ*norm(v.-x)^2<=c
                        break
                    end

                    αx*=ρ

                    if isnan(αx) || αx<10^-17
                        return x_, histF, histnψ, ls
                    end
                end
                if Fz>Fv
                    x=v
                    Fx=Fv
                    x_2=x
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
                return x_, histF, histnψ, ls
            end
        end
        t1=time()
        ∂fx=∂f(x)
        maybe=t1-time()

        elapsed=t1-start
        nψ=norm(∂fx.-∂fx_.+(x_2.-x).*αx_, p)
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

    return x, histF, histnψ, ls
end