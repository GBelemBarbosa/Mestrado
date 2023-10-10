using Plots
using LaTeXStrings

function mirror_C(f:: Function, ∂f:: Function, ∇ω:: Function, step:: Function, tₖ:: Function, dual_norm:: Function, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps) 
    x=x₀
    hist=[f(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)
        tₖ=tₖ(k, ∂fx)
        x=step(x, tₖ.*∂fx.-∇ω(x), tₖ)

        push!(hist, f(x))

        if dual_norm(∂fx)<ϵ
            break
        end
    end 

    println(dual_norm(∂f(x)), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end