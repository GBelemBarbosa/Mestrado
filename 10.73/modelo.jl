function Λ∂f(x:: Vector{<:Number}, ∂f:: Function)
    ∂fx=∂f(x)
    iₖ=argmax(∂fx)

    return sign(∂fx[iₖ]).*Float64.(I[1:n, iₖ])
end

function Lₖ(L:: Number, k:: Int64, x:: Vector{<:Number}, ∂fx⃰:: Vector{<:Number}, d∂x:: Float64; η=η, γ=γ, f=f, dual_norm=dual_norm) #Backtracking procedure B4
    fx=f(x)
    step=d∂x/L
    aux=(γ/L)*d∂x^2
    z=x.-step.*∂fx⃰

    while fx-f(z)<aux
        aux/=η
        step/=η
        L*=η
        z=x.-step.*∂fx⃰
    end

    return L, z
end

include("../Métodos/Non-Euclidean methods/Non-Euclidean_subgradient_descent_plot.jl")