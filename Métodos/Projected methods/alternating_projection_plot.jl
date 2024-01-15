using Plots
using LaTeXStrings

function alternating_projection(PS₁:: Function, PS₂:: Function, d₁:: Function, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps)  
    x=x₀
    hist=[d₁(x)]
    
    for k=0:k_max
        x=PS₂(PS₁(x)) #Projeta xₖ em S₁ e S₂ consecutivamente
        dx=d₁(x)

        push!(hist, dx)

        if dx<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"d_1(x^{(k)})",
                label=false)
end