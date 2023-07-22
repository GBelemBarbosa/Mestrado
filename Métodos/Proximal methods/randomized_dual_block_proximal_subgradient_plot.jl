using Plots
using LaTeXStrings

function RDBPG(f:: Function, G:: Function, step:: Function, Lₖ:: Function, y₀:: Array{Number, N}, s:: Vector{Number}, p:: Int64, n:: Int64, k_max:: Int64) where {N}
    y=y₀
    x=sumy=sum(y[i*n+1:(i+1)*n] for i=0:p-1)
    x_=x
    L=s
    hist=Float64[]
    
    for k=0:k_max
        iₖ=rand(1:p)
        x_, x=x, step(sumy)

        push!(hist, f(x)+G(x))

        L[iₖ], proxgᵢ=Lₖ(L[iₖ], iₖ, k, y, x) #Backtracking mais prox computation
        sumy.-=y[(iₖ-1)*n+1:i*n]
        aux=y[(iₖ-1)*n+1:i*n].+=(proxgᵢ.-x)./L[iₖ]
        sumy.+=aux
    end 

    println(norm(x.-x_, Inf), " ", f(x)+G(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end