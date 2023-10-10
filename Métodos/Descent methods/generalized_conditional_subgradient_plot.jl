using Plots
using LaTeXStrings

function generalized_conditional_subgradient(f:: Function, g:: Function, ∂f:: Function, ∂:: Function, tₖ:: Function, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps, q=Inf) 
    x=x₀
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)
        p=∂(x, ∂fx)     
        x.+=tₖ(k, x, p, ∂fx).*(p.-x)

        push!(hist, f(x)+g(x))

        if norm(p, q)<ϵ
            break
        end
    end 

    println(norm(∂(x, ∂f(x)), q), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end