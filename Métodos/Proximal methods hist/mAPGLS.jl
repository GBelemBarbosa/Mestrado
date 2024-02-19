using Plots
using LaTeXStrings

function mAPGLS(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, ρ:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    v=xsₖ=ysₖ=y=z=x_=x=x₀
    nysₖ=nxsₖ=αy=αx=t=1.0
    Fv=Fz=Fx=F(x)
    ∂fx=∂f(x)
    hist=[Fx]
    histnψ=[]
    ls=0
    
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
                return x_, hist, histnψ, ls
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
                return x_, hist, histnψ, ls
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
        ∂fx=∂f(x)

        push!(hist, Fx)
        nψ=norm(∂fx.-∂fx_.+(x_2.-x)./αx_, p)
        push!(nψ, histnψ)

        if nψ<ϵ || k==k_max
            break
        end
        k+=1

        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
        αy=nysₖ/(ysₖ'*(∂f(z).-∂fy))
        αx=nxsₖ/(xsₖ'*(∂f(v).-∂fx))
    end 

    return x, hist, histnψ, ls
end