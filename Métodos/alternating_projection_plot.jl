using Plots
using LaTeXStrings

function alternating_projection(PS₁:: Function, PS₂:: Function, d₁:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: K)  where {T, N, K}
    x=x₀
    hist=[d₁(x)]
    
    for k=0:k_max
        x=PS₂(PS₁(x)) #Projeta xₖ em S₁ e S₂ consecutivamente
        dx=d₁(x)

        push!(hist, dx)

        if dx<ϵ
            break
        end
    end 

    println(d₁(x))
    scatter(eachindex(hist), hist, 
                title=L"d_1(x^{(k)})",
                label=false)
end