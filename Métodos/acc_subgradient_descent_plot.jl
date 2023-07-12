using Plots
using LaTeXStrings

function acc_subgradient_descent(f:: Function, ∂f:: Function, β:: Number, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K) where {T, N, K}
    x₋, x=x₀, x₀
    a₋, a=1, 1
    u=0
    hist=[f(x)]
    
    for k=0:k_max        
        u=x.+a*(1/a₋-1).*(x.-x₋)
        ∂fu=∂f(u)

        if norm(∂fu, Inf)<ϵ
            break
        end

        x₋, x=x, u.-∂fu./β
        a₋, a=a, (sqrt(a^4+4*a^2)-a^2)/2

        push!(hist, f(x))
    end 

    println(norm(∂f(x), Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                yscale=:log10,
                title=L"f(x^{(k)})\ (accelerated)",
                label=false)
end