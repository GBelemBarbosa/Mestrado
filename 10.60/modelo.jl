using LinearAlgebra

include("variables.jl")

f(x:: Vector{<:Number}; A=A, b=b)=norm(A*x.-b, 2)^2/2

g(x:: Vector{<:Number}; λ=λ)=λ*norm(x, 1)

h(x:: Vector{<:Number}; D=D)=norm(D*x, 1)

function hμ(x:: Vector{<:Number}; D=D, μ=μ, p=p)
    sum=0.0
    y=D*x

    for i=1:p
        ny=abs(y[i])
        
        if ny>μ
            sum+=ny-μ/2
        else 
            sum+=(ny^2)/(2*μ)
        end
    end

    return sum
end

Τ(λ:: Number, x:: Vector{<:Number})=[max(abs(x[i])-λ, 0)*sign(x[i]) for i=1:length(x)] #Soft threshholding

function ∇Fμ(x:: Vector{<:Number}; A=A, b=b, D=D, μ=μ, Τ=Τ)
    Dx=D*x

    return A'*(A*x.-b).+D'*(Dx.-Τ(μ, Dx))./μ
end

proxgL(L, x:: Vector{<:Number}, λ=λ, Τ=Τ)=Τ(λ/L, x)
    
include("../Métodos/Proximal methods/S-FISTA_plot.jl")

x, p=SFISTA(f, g, h, hμ, ∇Fμ, proxgL, Lf, μ, α, x₀, k_max, eps())