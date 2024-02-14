using Plots
using LaTeXStrings

function mAPG(F:: Function, ∂f:: Function, proxα:: Function, x₀:: Array{<:Number}, αx:: Number, αy:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=z=x_=x=x₀
    t=1
    
    k=1
    while true
        ∂fy=∂f(y)
        z=proxα(αy, y, ∂fy)
        ∂fx=∂f(x)
        v=proxα(αx, x, ∂fx)
        if F(z)>F(v)
            x=v
        else
            x=z
        end

        if norm(∂fy, p)<ϵ || k==k_max
            break
        end
        k+=1

        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+(t_/t).*(z.-x).+((t_-1)/t).*(x.-x_)
        x_=x
    end 

    return x
end