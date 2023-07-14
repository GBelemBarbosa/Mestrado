using Plots
using LaTeXStrings

function incremental_projected_subgradient(f:: Function, ∂f:: Function, PC:: Function, m:: Int64, tₖ:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K) where {T, N, K}
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        ∂fmax=0.0

        for i=1:m
            ∂fxᵢ=∂f(x, i)
            aux=norm(∂fxᵢ, Inf)
            if ∂fmax<aux
                ∂fmax=aux
            end
            x=PC(x.-tₖ(k, ∂f).*∂fxᵢ)

            push!(hist, f(x))
        end

        if ∂fmax<ϵ
            break
        end
    end 

    println(max([norm(∂f(x, i), Inf) for i=1:m]), " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k, i)})",
                label=false)
end