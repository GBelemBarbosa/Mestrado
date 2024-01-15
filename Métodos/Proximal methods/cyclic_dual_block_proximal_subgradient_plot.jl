using Plots
using LaTeXStrings

function CDBPG(f:: Function, g:: Function, step:: Function, Lₖ:: Function, y₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, n:: Int64, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=y₀
    x=sum_y=sum(y[i*n+1:(i+1)*n] for i=0:p-1)
    x_=copy(x)
    L=s
    hist=[f(x)+g(x)]
    
    for k=0:k_max
        n∂x=0.0

        for i=1:p
            x_, x=x, step(sum_y)

            push!(hist, f(x)+g(x))

            aux=norm(x.-x_, p)            
            if n∂x<aux
                n∂x=aux
            end 
            L[i], proxgᵢ=Lₖ(L[i], i, k, y, x) #Backtracking mais prox computation
            sum_y.-=y[(i-1)*n+1:i*n]
            aux=y[(i-1)*n+1:i*n].+=(proxgᵢ.-x)./L[i]
            sum_y.+=aux
        end

        if n∂x<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(norm(x.-x_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end