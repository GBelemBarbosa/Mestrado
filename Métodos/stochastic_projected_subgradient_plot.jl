using Plots
using LaTeXStrings

function projected_subgradient(f:: Function, E∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{T, N}, k_max:: Int64) where {T, N}
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        x=PC(x.-tₖ(k).*E∂f(x))

        push!(hist, f(x))
    end 

    println(x, " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end