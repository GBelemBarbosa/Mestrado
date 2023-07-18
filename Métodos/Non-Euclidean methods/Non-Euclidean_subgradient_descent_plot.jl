using Plots
using LaTeXStrings

function NE_subgradient_descent(f:: Function, Λ∂f:: Function, Lₖ:: Function, dual_norm:: Function, x₀:: Array{T, N}, s:: Number, k_max:: Int64, ϵ:: Number) where {T, N}
    x=x₀
    L=s
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx⃰=Λ∂f(x)
        d∂x=dual_norm(∂fx⃰)

        if d∂x<ϵ
            break
        end

        L, x=Lₖ(L, k, x, ∂fx⃰, d∂x) #Backtracking mais atualização

        push!(hist, f(x))
    end 

    println(norm(Λ∂f(x), Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end