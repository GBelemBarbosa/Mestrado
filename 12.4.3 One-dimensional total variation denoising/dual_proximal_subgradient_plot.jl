using Plots
using LaTeXStrings

#Versão para operador multiplicação por matriz
function dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, Α:: Array{K, M}, Lₖ:: Function, y₀:: Array{T, N}, s:: Number, k_max:: Int64, ϵ:: Number) where  {K, M, T, N}
    y=y₀
    x=A'*y
    x_=x
    L=s
    hist=Float64[]
    
    for k=0:k_max
        x_, x=x, step(Α'*y)

        push!(hist, f(x)+g(A*x))
        
        if norm(x.-x_, Inf)<ϵ
            break
        end 

        Αx=Α*x
        L, prox=Lₖ(L, k, y, Αx) #Backtracking mais prox computation
        y=y.+(prox.-Αx)./L
    end 

    println(norm(x.-x_, Inf), " ", f(x)+g(A*x))
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end