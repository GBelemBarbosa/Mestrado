using Plots
using LaTeXStrings
using Random

function randomized_proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{T, N}, s:: Vector{Number}, p:: Int64, k_max:: Int64) where {T, N}
    x=x₀
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        iₖ=rand(1:p)

        ∂fxᵢ=∂f(x, iₖ)
        fx=f(x)
        L[iₖ], x=Lₖ(L[iₖ], iₖ, k, x, fx, ∂fxᵢ) #Backtracking mais atualização

        push!(hist, fx+g(x))
    end 

    println(max([norm(∂f(x, i), Inf) for i=1:p]), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end