using Plots
using LaTeXStrings

function PG(F:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, L₀:: Number, k_max:: Int64; ϵ=eps(), p=Inf)
    Fx=F(x₀)
    start=time()
    x_=x=x₀
    L=L₀
    ∂fx=∂f(x)
    histnψ=Tuple{Float64, Float64}[]
    histF=[(time()-start, Fx)]
    
    k=1
    while true
        L, x=Lₖ(L, k, x, ∂fx) #Backtracking mais atualização
        ∂fx_, ∂fx=∂fx, ∂f(x)

        t1=time()
        elapsed=t1-start
        nψ=norm(∂fx.-∂fx_.+(x_.-x).*L, p)
        push!(histF, (elapsed, F(x)))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1

        x_=x
    end 

    return x, histF, histnψ
end