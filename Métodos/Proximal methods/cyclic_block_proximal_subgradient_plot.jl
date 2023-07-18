using Plots
using LaTeXStrings

function proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{T, N}, s:: Vector{Number}, p:: Int64, k_max:: Int64, ϵ:: Number) where {T, N}
    x=x₀
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fmax=0.0

        for i=1:p
            ∂fxᵢ=∂f(x, i)
            aux=norm(∂fxᵢ, Inf)
            if ∂fmax<aux
                ∂fmax=aux
            end
            fx=f(x)
            L[i], x=Lₖ(L[i], i, k, x, fx, ∂fxᵢ) #Backtracking mais atualização

            push!(hist, fx+g(x))
        end

        if ∂fmax<ϵ
            break
        end
    end 

    println(max([norm(∂f(x, i), Inf) for i=1:p]), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end