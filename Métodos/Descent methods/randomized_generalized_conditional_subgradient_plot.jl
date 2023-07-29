using Plots
using LaTeXStrings

function RGBCG(f:: Function, g:: Function, ∂f:: Function, ∂:: Function, tₖ:: Function, x₀:: Array{<:Number}, p:: Int64, block_index:: Vector{UnitRange{Int64}}, k_max:: Int64) 
    x=x₀
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        iₖ=rand(1:p)
        ∂fxᵢ=∂f(x, iₖ)
        pᵢ=∂(x, ∂fxᵢ, iₖ)
        
        index=block_index[iₖ]
        x[index].+=tₖ(k, x, pᵢ, ∂fxᵢ).*(pᵢ.-x[index])

        push!(hist, f(x)+g(x))
    end 

    println(max([norm(∂(x, ∂f(x, i), i), Inf) for i=1:p]), " ", f(x)+g(x))
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end