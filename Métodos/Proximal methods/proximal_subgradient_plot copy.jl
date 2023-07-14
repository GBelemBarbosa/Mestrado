using Plots
using LaTeXStrings

function mirror_C(f:: Function, g:: Function, ∂f:: Function, prox_gLₖ:: Function, Lₖ:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K) where {T, N, K}
    x=x₀
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)

        if norm(∂fx, Inf)<ϵ
            break
        end

        Lₖ=Lₖ(Lₖ, k, x, ∂fx)
        x=prox_gLₖ(Lₖ, x.-∂fx./Lₖ)

        push!(hist, f(x)+g(x))
    end 

    println(dual_norm(∂f(x)), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end