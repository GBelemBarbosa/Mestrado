using Plots
using LaTeXStrings

function mAPG(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, αx:: Number, αy:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=z=x=x₀
    ∂fx=∂f(x)
    t=1.0
    hist=[F(x)]
    histnψ=[]
    
    k=1
    while true
        ∂fy=∂f(y)
        z=proxα(αy, y, ∂fy)
        v=proxα(αx, x, ∂fx)
        Fz=F(z)
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
        ∂fx=∂f(x)

        push!(hist, Fx)
        nψ=norm(∂fx.-∂fx_.+(x_2.-x).*αx_, p)
        push!(nψ, histnψ)

        if nψ<ϵ || k==k_max
            break
        end
        k+=1

        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
    end 

    return x, hist, histnψ
end