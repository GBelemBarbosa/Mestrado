using Plots
using LaTeXStrings

function VFISTA(f:: Function, g:: Function, ∇f:: Function, prox_gLf:: Function, Lf:: Number, σ:: Number, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=x=x₀
    rκ=sqrt(Lf/σ)
    step=(rκ-1)/(rκ+1)
    hist=[f(x)+g(x)]
    
    k=1
    while true
        ∇fy=∇f(y)
        x_, x=x, prox_gLf(y.-∇fy./Lf)
        push!(hist, f(x)+g(x))

        if norm(∇fy.-∇f(x).+(y.-x).*Lf, p)<ϵ || k==k_max
            break
        end
        k+=1

        y=x.+step.*(x.-x_)
    end 

    println(dual_norm(∇f(x)), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end