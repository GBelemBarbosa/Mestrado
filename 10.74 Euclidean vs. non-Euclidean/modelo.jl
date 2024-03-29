using LinearAlgebra

include("variables.jl")

f(x:: Vector{<:Number}; A=A, b=b)=x'*A*x/2+b'*x

∂f(x:: Vector{<:Number}; A=A, b=b)=A*x.+b

tₖ(k:: Int64, ∂f:: Vector{<:Number}; Lf₂=Lf₂)=1/Lf₂

include("../Métodos/Descent methods/subgradient_descent_plot.jl")

x, p₁=subgradient_descent(f, ∂f, tₖ, copy(x₀), Int64(ceil(k_max/n))) #Para entender o k_max/n, ver Remark 10.75

function Λ∂f(x:: Vector{<:Number}; A=A, b=b, n=n)
    maxabs, max=0.0, 0.0
    iₖ=0
    
    for i=1:n
        aux=A[i, :]'*x+b[i]
        absaux=abs(aux)
        if absaux>maxabs
            iₖ=i
            max=aux
            maxabs=absaux
        end
    end
    vec=zeros(Float64, n)
    vec[iₖ]=max

    return vec
end

Lₖ(L:: Number, k:: Int64, x:: Vector{<:Number}, ∂fx⃰:: Vector{<:Number}, d∂x:: Float64)=L, x.-∂fx⃰./L

dual_norm(x:: Vector{<:Number})=norm(x, Inf)

include("../Métodos/Non-Euclidean methods/Non-Euclidean_subgradient_descent_plot.jl")

x, p₂=NE_subgradient_descent(f, Λ∂f, Lₖ, dual_norm, x₀, Lf₁, k_max)

plot(p₁, p₂)