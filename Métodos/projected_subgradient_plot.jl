using Plots
using LaTeXStrings

function projected_subgradient(f:: Function, ∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K) where {T, N, K}
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)

        if norm(∂fx, Inf)<ϵ
            break
        end

        x=PC(x.-tₖ(k).*∂fx)

        push!(hist, f(x))
    end 

    println(norm(∂f(x), Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end