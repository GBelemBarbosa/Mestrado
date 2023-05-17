using Plots
using LaTeXStrings

function gradient_descent(f:: Function, ∇f:: Function, β:: Number, x₀, t_max:: Int64)
    x=x₀
    hist=[f(x)]
    
    for t=0:t_max
        ∇fx=∇f(x)
        if norm(∇fx, Inf)<ϵ
            break
        end
        x.-=∇fx./β

        push!(hist, f(x))
    end 

    println(norm(∇f(x), Inf), " ", f(x))
    scatter(eachindex(hist), hist, 
                yscale=:log10,
                title=L"f(x^{(k)})",
                label=false)
end