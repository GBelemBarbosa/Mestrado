using Plots
using LaTeXStrings

function stochastic_projected_subgradient(f:: Function, E∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K) where {T, N, K}
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        x_, x=x, PC(x.-tₖ(k, ∂f).*E∂f(x))

        push!(hist, f(x))

        if norm(x.-x_, Inf)<ϵ
            break
        end
    end 

    println(x, " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end