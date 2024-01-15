using Plots
using LaTeXStrings
using Random

function randomized_proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, k_max:: Int64; q=Inf) 
    x=x₀
    L=s
    fx=f(x)
    hist=[fx+g(x)]
    
    k=1
    while true
        iₖ=rand(1:p)

        ∂fxᵢ=∂f(x, iₖ)
        L[iₖ], x=Lₖ(L[iₖ], iₖ, k, x, fx, ∂fxᵢ) #Backtracking mais atualização
        fx=f(x)

        push!(hist, fx+g(x))

        if k==k_max
            break
        end
        k+=1
    end 

    println(max([norm(∂f(x, i), q) for i=1:p]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end