using Plots
using LaTeXStrings

function alternating_minimization(F:: Function, minimize:: Function, x₀:: Array{<:Number}, block_index:: Vector{UnitRange{Int64}}, k_max:: Int64; ϵ=eps, p=Inf) 
    x=x₀
    x_=copy(x)
    Fx=F(x)
    hist=[Fx]
    
    for k=0:k_max
        index=block_index[1]
        x_, x[index]=x, minimize(x, 1)
        n∂xᵢ=norm(x[index].-x_[index], p)
        Fx+=F(x, 1)-F(x_, 1)

        push!(hist, Fx)

        index=block_index[2]
        x_, x[index]=x, minimize(x, 2)
        aux=norm(x[index].-x_[index], p)
        if n∂xᵢ<aux
            n∂xᵢ=aux
        end
        Fx+=F(x, 2)-F(x_, 2)

        push!(hist, Fx)

        if n∂xᵢ<ϵ
            break
        end
    end 

    println(max([norm(x[i].-x_[i], p) for i=block_index]), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end