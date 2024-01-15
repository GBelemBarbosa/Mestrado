using Plots
using LaTeXStrings

function stochastic_projected_subgradient(f:: Function, E∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{<:Number}, k_max:: Int64) 
    x=x₀
    hist=[f(x)]
    
    k=1
    while true
        E∂fx=E∂f(x)
        x=PC(x.-tₖ(k, E∂fx).*E∂fx)

        push!(hist, f(x))

        if k==k_max
            break
        end
        k+=1
    end 

    println(x, " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end