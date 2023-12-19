using Plots
using LaTeXStrings

function projected_subgradient(f:: Function, ∂f:: Function, PC:: Function, tₖ:: Function, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps(), p=Inf) 
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)
        x=PC(x.-tₖ(k, ∂fx).*∂fx)

        push!(hist, f(x))

        if norm(∂fx, p)<ϵ
            break
        end
    end 

    println(norm(∂f(x), p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end