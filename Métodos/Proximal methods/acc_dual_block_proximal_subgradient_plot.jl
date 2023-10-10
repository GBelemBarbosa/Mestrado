using Plots
using LaTeXStrings

#Versão para operador multiplicação por matriz
function ADBPS(f̃:: Function, g₂:: Function, step:: Function, Lₖ:: Function, y₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps, p=Inf) 
    w=x=y=y₀
    u_=u=A'*y
    t=1
    L=s
    hist=[f̃(x)+g₂(x)]
    
    for k=0:k_max
        x=step(y)
        u_, u=u, step(w)

        push!(hist, f̃(x)+g₂(x))

        if norm(u.-u_, p)<ϵ
            break
        end 

        L, prox=Lₖ(L, k, w, u) #Backtracking + prox computation
        y_, y=y, w.+(prox.-u)./L
        t_, t=t, (1+sqrt(1+4*t^2))/2
        w=y.+((t_-1)/t).*(y.-y_) 
    end 

    println(norm(u.-u_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end