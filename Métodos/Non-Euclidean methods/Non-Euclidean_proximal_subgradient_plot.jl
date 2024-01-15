using Plots
using LaTeXStrings

function NE_proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, ∇ω:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    x=x₀
    L=s
    hist=[f(x)+g(x)]
    
    k=1
    while true
        ∂fx=∂f(x)
        ∇ωx=∇ω(x)
        fx=f(x)
        L, x=Lₖ(L, k, x, fx, ∂fx, ∇ωx) #Backtracking mais atualização

        push!(hist, fx+g(x))

        if norm(∂fx, p)<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(dual_norm(∂f(x)), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end