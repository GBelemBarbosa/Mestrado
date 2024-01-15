using Plots
using LaTeXStrings

function RGBCG(f:: Function, g:: Function, ∂f:: Function, ∂:: Function, tₖ:: Function, x₀:: Array{<:Number}, p:: Int64, block_index:: Vector{UnitRange{Int64}}, k_max:: Int64; p=Inf) 
    x=x₀
    hist=[f(x)+g(x)]
    
    k=1
    while true
        iₖ=rand(1:p)
        ∂fxᵢ=∂f(x, iₖ)
        pᵢ=∂(x, ∂fxᵢ, iₖ)
        
        index=block_index[iₖ]
        x[index].+=tₖ(k, x, pᵢ, ∂fxᵢ).*(pᵢ.-x[index])

        push!(hist, f(x)+g(x))

        if k==k_max
            break
        end
        k+=1
    end 

    println(max([norm(∂(x, ∂f(x, i), i), p) for i=1:p]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end