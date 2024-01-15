using Plots
using LaTeXStrings

function NE_subgradient_descent(f:: Function, Λ∂f:: Function, Lₖ:: Function, dual_norm:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    x=x₀
    L=s
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx⃰=Λ∂f(x)
        d∂x=dual_norm(∂fx⃰)
        L, x=Lₖ(L, k, x, ∂fx⃰, d∂x) #Backtracking mais atualização

        push!(hist, f(x))

        if d∂x<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(norm(Λ∂f(x), p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end