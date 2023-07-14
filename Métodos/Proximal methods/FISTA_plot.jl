using Plots
using LaTeXStrings

function mirror_C(f:: Function, g:: Function, ∂f:: Function, prox_gLₖ:: Function, Lₖ:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K) where {T, N, K}
    y, x=x₀, x₀
    t=1
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        ∂fx=∂f(x)

        if norm(∂fx, Inf)<ϵ
            break
        end

        Lₖ=Lₖ(Lₖ, k, y, ∂f(y))
        x_, x=x, prox_gLₖ(Lₖ, x.-∂fx./Lₖ)
        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/(t)).*(x.-x_)

        push!(hist, f(x)+g(x))
    end 

    println(dual_norm(∂f(x)), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end