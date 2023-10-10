using Plots
using LaTeXStrings

#Versão para operador multiplicação por matriz
function fast_dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, Α:: Array{<:Number, M}, Lₖ:: Function, y₀:: Array{<:Number, N}, s:: Number, k_max:: Int64; ϵ=eps, p=Inf) where {M, N}
    w=y=y₀
    u=A'*y
    u_=u
    t=1
    L=s
    hist=[f(x)+g(A*x)]
    
    for k=0:k_max
        u_, u=u, step(Α'*w)
        x=step(Α'*y)

        push!(hist, f(x)+g(A*x))
        
        if norm(u.-u_, p)<ϵ
            break
        end 

        Αu=Α*u
        L, prox=Lₖ(L, k, w, Αu) #Backtracking mais prox computation
        y_, y=y, w.+(prox.-Αu)./L
        t_, t=t, (1+sqrt(1+4*t^2))/2
        w=y.+((t_-1)/t).*(y.-y_)
    end 

    println(norm(u.-u_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end

#Versão operador genérico
function fast_dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, Α:: Function, ΑT:: Function, Lₖ:: Function, y₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps, p=Inf) 
    w, y=y₀, y₀
    u=ΑT(y)
    u_=u
    t=1
    L=s
    hist=[f(x)+g(Α(x))]
    
    for k=0:k_max
        u_, u=u, step(ΑT(w))
        x=step(ΑT(y))

        push!(hist, f(x)+g(Α(x)))
        
        if norm(u.-u_, p)<ϵ
            break
        end 

        Αu=Α(u)
        L, prox=Lₖ(L, k, w, Αu) #Backtracking mais prox computation
        y_, y=y, w.+(prox.-Αu)./L
        t_, t=t, (1+sqrt(1+4*t^2))/2
        w=y.+((t_-1)/t).*(y.-y_)
    end 

    println(norm(u.-u_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end