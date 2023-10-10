using Plots
using LaTeXStrings

function acc_subgradient_descent(f:: Function, ∂f:: Function, β:: Number, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps, p=Inf) 
    x₋=x=x₀
    a₋=a=1
    u=0
    hist=[f(x)]
    
    for k=0:k_max        
        u=x.+a*(1/a₋-1).*(x.-x₋)
        ∂fu=∂f(u)
        x₋, x=x, u.-∂fu./β
        a₋, a=a, (sqrt(a^4+4*a^2)-a^2)/2

        push!(hist, f(x))

        if norm(∂fu, p)<ϵ
            break
        end
    end 

    println(norm(∂f(x), p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})\ (accelerated)",
                label=false)
end