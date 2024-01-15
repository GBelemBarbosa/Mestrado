using Plots
using LaTeXStrings

#Versão para operador multiplicação por matriz
function dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, Α:: Array{<:Number, M}, Lₖ:: Function, y₀:: Array{<:Number, N}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) where {M, N}
    y=y₀
    x=Α'*y
    x_=x
    L=s
    hist=[f(x)+g(A*x)]
    
    for k=0:k_max
        x_, x=x, step(Α'*y)

        push!(hist, f(x)+g(A*x))
        
        if norm(x.-x_, p)<ϵ || k==k_max
            break
        end
        k+=1 

        Αx=Α*x
        L, prox=Lₖ(L, k, y, Αx) #Backtracking mais prox computation
        y.+=(prox.-Αx)./L
    end 

    println(norm(x.-x_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end

#Versão operador genérico
function dual_proximal_subgradient(f:: Function, g:: Function, step:: Function, Α:: Function, ΑT:: Function, Lₖ:: Function, y₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=y₀
    x=ΑT(y)
    x_=x
    L=s
    hist=[f(x)+g(Α(x))]
    
    for k=0:k_max
        x_, x=x, step(ΑT(y))

        push!(hist, f(x)+g(Α(x)))
        
        if norm(x.-x_, p)<ϵ || k==k_max
            break
        end
        k+=1 

        Αx=Α(x)
        L, prox=Lₖ(L, k, y, Αx) #Backtracking mais prox computation
        y.+=(prox.-Αx)./L
    end 

    println(norm(x.-x_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end