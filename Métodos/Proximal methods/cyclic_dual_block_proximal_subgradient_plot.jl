using Plots
using LaTeXStrings

function CDBPG(f:: Function, g:: Function, step:: Function, Lₖ:: Function, y₀:: Array{<:Number}, s:: Vector{<:Number}, p:: Int64, n:: Int64, k_max:: Int64, ϵ:: Number) where {N}
    y=y₀
    x=sumy=sum(y[i*n+1:(i+1)*n] for i=0:p-1)
    x_=x
    L=s
    hist=Float64[]
    
    for k=0:k_max
        ∂xmax=0.0

        for i=1:p
            x_, x=x, step(sumy)

            push!(hist, f(x)+g(x))

            aux=norm(x.-x_, Inf)            
            if ∂xmax<aux
                ∂xmax=aux
            end 
            L[i], proxgᵢ=Lₖ(L[i], i, k, y, x) #Backtracking mais prox computation
            sumy.-=y[(i-1)*n+1:i*n]
            aux=y[(i-1)*n+1:i*n].+=(proxgᵢ.-x)./L[i]
            sumy.+=aux
        end

        if ∂xmax<ϵ
            break
        end
    end 

    println(norm(x.-x_, Inf), " ", f(x)+g(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end