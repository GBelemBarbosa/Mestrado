using Plots
using LaTeXStrings

function nmAPG(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, αx:: Number, αy:: Number, η:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=z=x_=x=x₀
    c=F(x)
    ∂fx=∂f(x)
    q=t=1.0
    
    k=1
    while true
        ∂fy=∂f(y)
        z=proxα(αy, y, ∂fy)
        Fz=F(z)
        if Fz<=c-δ*norm(z.-y)^2
            x=z
            Fx=Fz
            x_2=y
            ∂fx_=∂fy
            αx_=αy
        else
            v=proxα(αx, x, ∂fx)
            Fv=F(v)
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
        end
        ∂fx=∂f(x)

        if norm(∂fx.-∂fx_.+(x_2.-x)./αx_, p)<ϵ || k==k_max
            break
        end
        k+=1
        
        t_, t=t, (1+sqrt(1+4*t^2))/2
        q_, q=q, η*q+1
        c=(η*q_*c+Fx)/q
        y=x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
    end 

    return x
end