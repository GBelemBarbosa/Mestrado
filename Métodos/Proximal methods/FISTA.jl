using Plots
using LaTeXStrings

function FISTA(∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=x=x₀
    t=1
    L=s
    
    k=1
    while true
        ∂fy=∂f(y)
        x_=x
        L, x=Lₖ(L, k, y, ∂fy) #Backtracking mais atualização
        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_)

        if norm(∂fy, p)<ϵ || k==k_max
            break
        end
        k+=1
    end 

    return x
end