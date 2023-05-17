using Plots
using LaTeXStrings

function acc_gradient_descent(f:: Function, ∇f:: Function, β:: Number, x₀, t_max:: Int64)
    x₋, x=x₀, x₀
    a₋, a=1, 1
    u=0
    hist=[f(x)]
    
    for t=0:t_max        
        u=x.+a*(1/a₋-1).*(x.-x₋)
        ∇fu=∇f(u)
        if norm(∇fu, Inf)<ϵ
            break
        end
        x₋, x=x, u.-∇fu./β

        push!(hist, f(x))

        a₋, a=a, (sqrt(a^4+4*a^2)-a^2)/2
    end 

    println(norm(∇f(u), Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                yscale=:log10,
                title=L"f(x^{(k)})\ (accelerated)",
                label=false)
end