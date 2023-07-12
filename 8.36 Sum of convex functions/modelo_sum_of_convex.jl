using LinearAlgebra

include("variables_sum.jl")

function f(x:: Vector{Float64}, i:: Int64)
    if i==1
        fx=sum(x)
    elseif i==2
        fx=norm(x, 2)
    else
        fx=norm(x, 1)
    end

    return fx
end

function Σfᵢ(x:: Vector{Float64}; m=m, f=f)
    sum=0.0

    for i=1:m
        sum+=f(x, i)
    end

    return sum
end

function ∂f(x:: Vector{Float64}, i:: Int64; n=n)
    if i==1
        fx=ones(n)
    elseif i==2
        nx=norm(x, 2)

        if nx!=0
            fx=x./nx
        else
            fx=randn(Float64, n)
            fx./=norm(fx, 2)
        end
    else
        fx=sign.(x)
    end
    
    return fx
end

function Σ∂fᵢ(x:: Vector{Float64}; m=m, n=n, ∂f=∂f) #Função objetivo
    sum=zeros(Float64, n)

    for i=1:m
        sum.+=∂f(x, i)
    end

    return sum
end

function E∂f(x:: Vector{Float64}; m=m, Θ=Θ, Σ∂fᵢ=Σ∂fᵢ)
    aux=Σ∂fᵢ(x)

    return (sqrt(2*Θ)/norm(aux)).*aux
end

PC(x:: Vector{Float64})=min.(max.(x, -1), 1) #Projeção de x em C

tₖ(k:: Int64)=1/sqrt(k+1) #Stepsize rule

include("../Métodos/stochastic_projected_subgradient_plot.jl")

p₁=projected_subgradient(Σfᵢ, E∂f, PC, tₖ, x₀, k_max)

E∂f(x:: Vector{Float64}; m=m, Θ=Θ, Lf=Lf, ∂f=∂f)=(sqrt(2*Θ)*m/Lf).*∂f(x, rand(1:m))

p₂=projected_subgradient(Σfᵢ, E∂f, PC, tₖ, x₀, k_max)

plot(p₁, p₂)