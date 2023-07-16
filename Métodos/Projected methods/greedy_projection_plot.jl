using Plots
using LaTeXStrings

function greedy_projection(PS:: Function, d:: Function, x₀:: Array{T, N}, k_max:: Int64, ϵ:: Number) where {T, N}
    x=x₀
    hist=Float64[]
    
    for k=0:k_max
        d_max, iₖ=d(x)

        push!(hist, d_max)

        if d_max<ϵ
            break
        end

        x=PS(x, iₖ) #Projeta xₖ em Sᵢₖ, sendo iₖ∈argmax(d(xₖ, Sᵢ))
    end 

    println(d(x))
    scatter(eachindex(hist), hist, 
                title=L"d_{max}(x^{(k)})",
                label=false)
end