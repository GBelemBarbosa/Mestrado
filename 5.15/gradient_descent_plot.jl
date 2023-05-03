using Plots
using LaTeXStrings

function gradient_descent(f:: Function, ∇f:: Function, β:: Float64, x₀, t_max:: Int64)
    x=x₀
    hist=Float64[]
    for t=1:t_max
        push!(hist, f(x))

        ∇fx=∇f(x)
        x-=∇fx/β

        if norm(∇fx, Inf)<ϵ
            break
        end
    end 

    println(norm(∇f(x), Inf))
    plot=scatter(eachindex(hist), hist, 
                yscale=:log10,
                title=L"Semilog-y plot of $f(x^{(k)})$",
                label=false)
end