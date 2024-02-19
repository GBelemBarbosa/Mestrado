using Plots
using LaTeXStrings

function CBPG(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, block_index:: Vector{UnitRange{Int64}}, k_max:: Int64; ϵ=eps(), q=Inf) 
    x=x₀
    L=s
    fx=f(x)
    hist=[fx+g(x)]
    
    k=1
    while true
        n∂x=0.0

        for i=1:p    
            x_=copy(x)        
            L[i], x[block_index[i]]=Lₖ(L[i], i, k, x, fx, ∂f(x, i)) #Backtracking mais atualização
            fx=f(x)

            aux=norm(x[block_index[i]].-x_[block_index[i]], q)
            if n∂x<aux
                n∂x=aux
            end

            push!(hist, fx+g(x))
        end

        if n∂x<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(max([norm(∂f(x, i), q) for i=1:p]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end