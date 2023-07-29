using Plots
using LaTeXStrings

function mirror_descent(f:: Function, ∂f:: Function, ∇ω:: Function, step:: Function, tₖ:: Function, dual_norm:: Function, x₀:: Array{<:Number}, k_max:: Int64, ϵ:: Number) 
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)

        if dual_norm(∂fx)<ϵ
            break
        end

        x=step(x, tₖ(k, ∂fx).*∂fx.-∇ω(x))

        push!(hist, f(x))
    end 

    println(dual_norm(∂f(x)), " ", f(x))
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end