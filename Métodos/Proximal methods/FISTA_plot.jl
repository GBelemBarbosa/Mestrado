using Plots
using LaTeXStrings

function FISTA(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64, ϵ:: Number) 
    y, x=x₀, x₀
    t=1
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fy=∂f(y)

        if norm(∂fy, Inf)<ϵ
            break
        end

        x_=x
        L, x=Lₖ(L, k, y, ∂fy) #Backtracking mais atualização
        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_)

        push!(hist, f(x)+g(x))
    end 

    println(norm(∂f(x), Inf), " ", f(x)+g(x))
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end