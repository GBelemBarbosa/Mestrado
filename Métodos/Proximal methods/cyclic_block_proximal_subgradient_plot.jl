using Plots
using LaTeXStrings

function CBPG(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, block_index:: Vector{UnitRange{Int64}}, k_max:: Int64; ϵ=eps(), q=Inf) 
    x=x₀
    L=s
    fx=f(x)
    hist=[fx+g(x)]
    
    for k=0:k_max
        n∂f=0.0

        for i=1:p
            ∂fxᵢ=∂f(x, i)
            aux=norm(∂fxᵢ, q)
            if n∂f<aux
                n∂f=aux
            end
            L[i], x[block_index[i]]=Lₖ(L[i], i, k, x, fx, ∂fxᵢ) #Backtracking mais atualização
            fx=f(x)

            push!(hist, fx+g(x))
        end

        if n∂f<ϵ
            break
        end
    end 

    println(max([norm(∂f(x, i), q) for i=1:p]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end