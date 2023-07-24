using Plots
using LaTeXStrings

function proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64, ϵ:: Number) where {N}
    x=x₀
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)

        if norm(∂fx, Inf)<ϵ
            break
        end

        fx=f(x)
        L, x=Lₖ(L, k, x, fx, ∂fx) #Backtracking mais atualização

        push!(hist, fx+g(x))
    end 

    println(norm(∂f(x), Inf), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end