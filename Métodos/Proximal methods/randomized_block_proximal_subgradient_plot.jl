using Plots
using LaTeXStrings
using Random

function randomized_proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, k_max:: Int64) 
    x=x₀
    L=s
    fx=f(x)
    hist=[fx+g(x)]
    
    for k=0:k_max
        iₖ=rand(1:p)

        ∂fxᵢ=∂f(x, iₖ)
        L[iₖ], x=Lₖ(L[iₖ], iₖ, k, x, fx, ∂fxᵢ) #Backtracking mais atualização
        fx=f(x)

        push!(hist, fx+g(x))
    end 

    println(max([norm(∂f(x, i), Inf) for i=1:p]), " ", f(x)+g(x))
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end