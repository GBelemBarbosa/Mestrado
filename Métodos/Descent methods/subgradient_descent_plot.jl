using Plots
using LaTeXStrings

function subgradient_descent(f:: Function, ∂f:: Function, tₖ:: Function, x₀:: Array{<:Number}, k_max:: Int64, ϵ:: Number)
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)

        if norm(∂fx, Inf)<ϵ
            break
        end
        
        x.-=tₖ(k, ∂fx).*∂fx

        push!(hist, f(x))
    end 

    println(norm(∂f(x), Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end