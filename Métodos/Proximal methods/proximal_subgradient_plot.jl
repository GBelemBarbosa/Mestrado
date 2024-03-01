using Plots
using LaTeXStrings

function proximal_subgradient(F:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, L₀:: Number, k_max:: Int64; ϵ=eps(), p=Inf)
    x=x₀
    L=L₀
    ∂fx=∂f(x)
    hist=[F(x)]
    
    k=1
    while true
        L, x=Lₖ(L, k, x, ∂fx) #Backtracking mais atualização
        ∂fx_, ∂fx=∂fx, ∂f(x)

        push!(hist, F(x))

        if norm(∂fx.-∂fx_.+(x_.-x).*L, p)<ϵ || k==k_max
            break
        end
        k+=1

        x_=x
    end 

    println(norm(∂f(x), p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end