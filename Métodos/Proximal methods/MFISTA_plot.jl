using Plots
using LaTeXStrings

function MFISTA(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=x=x₀
    t=1
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fy=∂f(y)
        x_=x
        L, z=Lₖ(L, k, y, ∂fy) #Backtracking mais atualização
        if f(x_)+g(x_)<f(z)+g(z)
            x=x_
        else
            x=z
        end

        push!(hist, f(x)+g(x))

        if norm(∂fy, p)<ϵ || k==k_max
            break
        end
        k+=1

        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_).+(t_/t).*(z.-x)
    end 

    println(dual_norm(∂f(x)), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end