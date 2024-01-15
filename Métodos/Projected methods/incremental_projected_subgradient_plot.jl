using Plots
using LaTeXStrings

function incremental_projected_subgradient(f:: Function, ∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{<:Number},  m:: Int64, k_max:: Int64; ϵ=eps(), q=Inf) 
    x=x₀
    hist=[f(x)]
    
    k=1
    while true
        n∂f=0.0

        for i=1:m
            ∂fxᵢ=∂f(x, i)
            aux=norm(∂fxᵢ, q)
            if n∂f<aux
                n∂f=aux
            end
            x=PC(x.-tₖ(k, ∂fxᵢ).*∂fxᵢ)

            push!(hist, f(x))
        end

        if n∂f<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(max([norm(∂f(x, i), q) for i=1:m]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k, i)})",
                label=false)
end