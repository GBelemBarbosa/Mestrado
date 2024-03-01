using Plots
using LaTeXStrings

function mAPGLS(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, ρ:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    Fx=F(x₀)
    start=time()
    v=xsₖ=ysₖ=y=z=x_=x=x₀
    nysₖ=nxsₖ=αy=αx=t=1.0
    Fv=Fz=Fx=F(x)
    ∂fx=∂f(x)
    histnψ=Tuple{Float64, Float64}[]
    ls=0
    histF=[(time()-start, Fx)]
    
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

            if Fz+δ*nysₖ<=Fy
                break
            end

            αy*=ρ

            if isnan(αy) || αy<10^-17
                return x_, histF, histnψ, ls
            end
        end

        while true
            v=proxα(αx, x, ∂fx)

            ls+=1

            Fv=F(v)
            xsₖ=v.-x
            nxsₖ=xsₖ'xsₖ

            if Fv+δ*nxsₖ<=Fx
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
            ∂fx_2=∂fx
            αx_=αx
        else
            x=z
            Fx=Fz
            x_2=y
            ∂fx_2=∂fy
            αx_=αy
        end
        ∂fx_, ∂fx=∂fx, ∂f(x)

        t1=time()
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
        y=x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
        αy=nysₖ/(ysₖ'*(∂f(z).-∂fy))
        αx=nxsₖ/(xsₖ'*(∂f(v).-∂fx_))
    end 

    return x, histF, histnψ, ls
end