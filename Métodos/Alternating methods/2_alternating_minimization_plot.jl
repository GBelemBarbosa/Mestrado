using Plots
using LaTeXStrings

function alternating_minimization(F:: Function, minimize:: Function, x₀:: Array{<:Number, N}, block_index:: Vector{<:Numberuple{Int64, Int64}}, k_max:: Int64, ϵ:: Number) where {N}
    x=x₀
    Fx=F(x)
    hist=[Fx]
    
    for k=0:k_max
        index=block_index[1][1]:block_index[1][2]
        x_, x[index]=x, minimize(x, 1)
        ∂xᵢmax=norm(x[index].-x_[index], Inf)
        Fx+=F(x, 1)-F(x_, 1)

        push!(hist, Fx)

        index=block_index[2][1]:block_index[2][2]
        x_, x[index]=x, minimize(x, 2)
        aux=norm(x[index].-x_[index], Inf)
        if ∂xᵢmax<aux
            ∂xᵢmax=aux
        end
        Fx+=F(x, 2)-F(x_, 2)

        push!(hist, Fx)

        if ∂xᵢmax<ϵ
            break
        end
    end 

    println(∂xᵢmax, " ", F(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end