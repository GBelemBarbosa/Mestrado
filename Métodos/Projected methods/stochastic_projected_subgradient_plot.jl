using Plots
using LaTeXStrings

function stochastic_projected_subgradient(f:: Function, E∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{<:Number, N}, k_max:: Int64) where {N}
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        E∂fx=E∂f(x)
        x=PC(x.-tₖ(k, E∂fx).*E∂fx)

        push!(hist, f(x))
    end 

    println(x, " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end