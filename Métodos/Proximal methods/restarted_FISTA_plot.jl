using Plots
using LaTeXStrings

include("FISTA.jl")

function restarted_FISTA(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64, N:: Int64, ϵ:: Number) where {N}
    z=x₀
    hist=[f(z)+g(z)]
    
    for k=0:k_max
        z=FISTA(∂f, Lₖ, z, s, N, ϵ)

        push!(hist, f(z)+g(z))

        if norm(∂f(z), Inf)<ϵ
            break
        end
    end 

    println(dual_norm(∂f(z)), " ", f(z)+g(z))
    scatter(eachindex(hist), hist, 
                title=L"F(x^{(N*k)})",
                label=false)
end