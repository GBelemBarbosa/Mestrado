using Plots
using LaTeXStrings

function fast_dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, Α:: Function, Lₖ:: Function, y₀:: Array{T, N}, s:: Number, k_max:: Int64, ϵ:: Number) where {T, N}
    w, y=y₀, y₀
    u=A'*y
    t=1
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        u_, u=u, step(Α'*w)
        x=step(Α'*y)

        push!(hist, f(x)+g(x))
        
        if norm(u.-u_, Inf)<ϵ
            break
        end 

        Αu=Α*u
        L, prox=Lₖ(L, k, w, Αu) #Backtracking mais prox computation
        y_, y=y, w.-Αu./L.+prox
        t_, t=t, (1+sqrt(1+4*t^2))/2
        w=y.+((t_-1)/t).*(y.-y_)
    end 

    x=step(Α'*y)
    println(norm(u.-u_, Inf), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end