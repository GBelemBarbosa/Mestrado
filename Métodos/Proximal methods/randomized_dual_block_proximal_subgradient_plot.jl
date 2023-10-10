using Plots
using LaTeXStrings

function RDBPG(f:: Function, G:: Function, step:: Function, Lₖ:: Function, y₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, n:: Int64, k_max:: Int64; q=Inf) 
    y=y₀
    x=sum_y=sum(y[i*n+1:(i+1)*n] for i=0:p-1)
    x_=copy(x)
    L=s
    hist=[f(x)+G(x)]
    
    for k=0:k_max
        iₖ=rand(1:p)
        x_, x=x, step(sum_y)

        push!(hist, f(x)+G(x))

        L[iₖ], proxgᵢ=Lₖ(L[iₖ], iₖ, k, y, x) #Backtracking mais prox computation
        sum_y.-=y[(iₖ-1)*n+1:i*n]
        aux=y[(iₖ-1)*n+1:i*n].+=(proxgᵢ.-x)./L[iₖ]
        sum_y.+=aux
    end 

    println(norm(x.-x_, q), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end