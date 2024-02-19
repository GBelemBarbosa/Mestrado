using Plots
using LaTeXStrings

function PG(F:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, L₀:: Number, k_max:: Int64; ϵ=eps(), p=Inf)
    x=x₀
    L=L₀
    ∂fx=∂f(x)
    hist=[F(x)]
    histnψ=[]
    
    k=1
    while true
        L, x=Lₖ(L, k, x, ∂fx) #Backtracking mais atualização
        ∂fx_, ∂fx=∂fx, ∂f(x)
        
        push!(hist, F(x))
        nψ=norm(∂fx.-∂fx_.+(x_.-x).*L, p)
        push!(nψ, histnψ)

        if nψ<ϵ || k==k_max
            break
        end
        k+=1
    end 

    return x, hist, histnψ
end