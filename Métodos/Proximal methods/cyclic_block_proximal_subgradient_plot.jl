using Plots
using LaTeXStrings

function CBPG(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{Number, N}, s:: Vector{Number}, p:: Int64, block_index:: Vector{Tuple{Int64, Int64}}, k_max:: Int64, ϵ:: Number) where {N}
    x=x₀
    L=s
    fx=f(x)
    hist=[fx+g(x)]
    
    for k=0:k_max
        ∂fmax=0.0

        for i=1:p
            ∂fxᵢ=∂f(x, i)
            aux=norm(∂fxᵢ, Inf)
            if ∂fmax<aux
                ∂fmax=aux
            end
            L[i], x[block_index[i][1]:block_index[i][2]]=Lₖ(L[i], i, k, x, fx, ∂fxᵢ) #Backtracking mais atualização
            fx=f(x)

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