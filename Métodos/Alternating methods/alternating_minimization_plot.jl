using Plots
using LaTeXStrings

function alternating_minimization(F:: Function, minimize:: Function, x₀:: Array{<:Number}, p:: Int64, block_index:: Vector{UnitRange{Int64}}, k_max:: Int64; ϵ=eps(), q=Inf) 
    x=x₀
    x_=copy(x)
    Fx=F(x)
    hist=[Fx]
    
    k=1
    while true
        n∂xᵢ=0.0

        for i=1:p
            index=block_index[i]
            x_, x[index]=x, minimize(x, i)
            aux=norm(x[index].-x_[index], q)
            if n∂xᵢ<aux
                n∂xᵢ=aux
            end
            Fx+=F(x[index], i)-F(x_[index], i)

            push!(hist, Fx)
        end

        if n∂xᵢif norm(∂fx, p)<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(max([norm(x[i].-x_[i], q) for i=block_index]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end