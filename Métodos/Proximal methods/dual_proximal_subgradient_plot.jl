using Plots
using LaTeXStrings

function dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, A:: Function, Lₖ:: Function, y₀:: Array{T, N}, s:: Number, k_max:: Int64, ϵ:: Number) where {T, N}
    y=y₀
    x=A'*y
    L=s
    hist=Float64[]
    
    for k=0:k_max
        x_, x=x, step(Α'*y)

        push!(hist, f(x)+g(x))
        
        if norm(x.-x_, Inf)<ϵ
            break
        end 

        Αx=Α*x
        L, prox=Lₖ(L, k, y, Αx) #Backtracking mais prox computation
        y=y.-Αx./L.+prox
    end 

    println(norm(x.-x_, Inf), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end