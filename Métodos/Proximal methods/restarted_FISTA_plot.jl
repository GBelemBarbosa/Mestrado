using Plots
using LaTeXStrings

include("FISTA.jl")

function restarted_FISTA(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64, N:: Int64; ϵ=eps(), p=Inf) 
    z=x₀
    hist=[f(z)+g(z)]
    
    k=1
    while true
        z=FISTA(∂f, Lₖ, z, s, N)

        push!(hist, f(z)+g(z))

        if norm(∂f(z), p)<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(norm(∂f(z), p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(N*k)})",
                label=false)
end