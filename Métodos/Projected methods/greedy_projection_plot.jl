using Plots
using LaTeXStrings

function greedy_projection(PS:: Function, d:: Function, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps) 
    x=x₀
    d_max, iₖ=d(x)
    hist=Float64[d_max]
    
    k=1
    while true
        x=PS(x, iₖ) #Projeta xₖ em Sᵢₖ, sendo iₖ∈argmax(d(xₖ, Sᵢ))
        
        push!(hist, d_max)

        if d_max<ϵ || k==k_max
            break
        end
        k+=1

        d_max, iₖ=d(x)
    end 

    println(hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"d_{max}(x^{(k)})",
                label=false)
end