using Plots
using LaTeXStrings

#Versão para operador multiplicação por matriz
function fast_dual_proximal_subgradient(f̃:: Function, g₂:: Function, step:: Function, Lₖ:: Function, y₀:: Array{<:Number, N}, s:: Number, k_max:: Int64, ϵ:: Number) where {N}
    w, y=y₀, y₀
    u=A'*y
    u_=u
    t=1
    L=s
    hist=Float64[]
    
    for k=0:k_max
        u_, u=u, step(w)
        x=step(y)

        push!(hist, f̃(x)+g₂(x))
        
        if norm(u.-u_, Inf)<ϵ
            break
        end 

        L, prox=Lₖ(L, k, w, u) #Backtracking mais prox computation
        y_, y=y, w.+(prox.-u)./L
        t_, t=t, (1+sqrt(1+4*t^2))/2
        w=y.+((t_-1)/t).*(y.-y_)
    end 

    x=step(y)
    println(norm(u.-u_, Inf), " ", f̃(x)+g₂(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end