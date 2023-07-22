using Plots
using LaTeXStrings

function alternating_minimization(F:: Function, minimize:: Function, x₀:: Array{Number, N}, p:: Int64, block_index:: Vector{Tuple{Int64, Int64}}, k_max:: Int64, ϵ:: Number) where {N}
    x=x₀
    Fx=F(x)
    hist=[Fx]
    
    for k=0:k_max
        ∂xᵢmax=0.0

        for i=1:p
            index=block_index[i][1]:block_index[i][2]
            x_, x[index]=x, minimize(x, i)
            aux=norm(x[index].-x_[index], Inf)
            if ∂xᵢmax<aux
                ∂xᵢmax=aux
            end
            Fx+=F(x[index], i)-F(x_[index], i)

            push!(hist, Fx)
        end

        if ∂xᵢmax<ϵ
            break
        end
    end 

    println(∂xᵢmax, " ", F(x))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(k, i)})",
                label=false)
end