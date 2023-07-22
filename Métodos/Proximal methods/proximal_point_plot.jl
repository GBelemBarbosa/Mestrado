using Plots
using LaTeXStrings

function proximal_point(g:: Function, prox_cg:: Function, c:: Number, x₀:: Array{Number, N}, k_max:: Int64, ϵ:: Number) where {N}
    x=x₀
    hist=[g(x)]
    
    for k=0:k_max
        x_, x=x, prox_cg(c, x)

        push!(hist, g(x))

        if norm(x.-x_, Inf)<ϵ
            break
        end
    end 

    println(norm(x.-x_, Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                title=L"g(x^{(k)})",
                label=false)
end