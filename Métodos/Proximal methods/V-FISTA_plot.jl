using Plots
using LaTeXStrings

function VFISTA(f:: Function, g:: Function, ∇f:: Function, prox_gLf:: Function, Lf:: Number, σ:: Number, x₀:: Array{T, N}, k_max:: Int64, ϵ:: Number) where {T, N}
    y, x=x₀, x₀
    rκ=sqrt(Lf/σ)
    step=(rκ-1)/(rκ+1)
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∇fx=∇f(x)

        if norm(∇fx, Inf)<ϵ
            break
        end

        x_, x=x, prox_gLf(x.-∇fx./Lf)
        y=x.+step.*(x.-x_)

        push!(hist, f(x)+g(x))
    end 

    println(dual_norm(∇f(x)), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end