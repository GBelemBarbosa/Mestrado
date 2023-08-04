include("variables_Ridge.jl")
using BenchmarkTools

f(x:: Vector{<:Number})=((A*x.-y)⋅(A*x.-y)+λ*x⋅x)/2
f(x:: Vector{<:Number}, A:: Matrix{Float64}, y:: Vector{<:Number}, λ:: Number)=((A*x.-y)⋅(A*x.-y)+λ*x⋅x)/2
fᵣ(x:: Vector{<:Number}; A=A, y=y, λ=λ)=((A*x.-y)⋅(A*x.-y)+λ*x⋅x)/2

g(x:: Vector{<:Number})=(norm(A*x.-y)^2+λ*norm(x)^2)/2
g(x:: Vector{<:Number}, A:: Matrix{Float64}, y:: Vector{<:Number}, λ:: Number)=(norm(A*x.-y)^2+λ*norm(x)^2)/2
gᵣ(x:: Vector{<:Number}; A=A, y=y, λ=λ)=(norm(A*x.-y)^2+λ*norm(x)^2)/2

function h(x:: Vector{<:Number})
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end
function h(x:: Vector{<:Number}, A:: Matrix{Float64}, y:: Vector{<:Number}, λ:: Number)
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end
function hᵣ(x:: Vector{<:Number}; A=A, y=y, λ=λ)
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end

f₂(x)=((A*x.-y)⋅(A*x.-y)+λ*x⋅x)/2
f₂(x, A, y, λ)=((A*x.-y)⋅(A*x.-y)+λ*x⋅x)/2

g₂(x)=(norm(A*x.-y)^2+λ*norm(x)^2)/2
g₂(x, A, y, λ)=(norm(A*x.-y)^2+λ*norm(x)^2)/2

function h₂(x)
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end
function h₂(x, A, y, λ)
    aux=A*x.-y
    return (aux⋅aux+λ*x⋅x)/2
end

function repeval(f:: Function)
    for i in 1:10000
        res=f(randn(n))
    end
end

function repeval(f:: Function, A:: Matrix{Float64}, y:: Vector{<:Number}, λ:: Number)
    for i in 1:10000
        res=f(randn(n), A, y, λ)
    end
end

@btime repeval(f)
@btime repeval(f, A, y, λ)
@btime repeval(fᵣ)
@btime repeval(f₂)
@btime repeval(f₂, A, y, λ)

@btime repeval(g)
@btime repeval(g, A, y, λ)
@btime repeval(gᵣ)
@btime repeval(g₂)
@btime repeval(g₂, A, y, λ)

@btime repeval(h)
@btime repeval(h, A, y, λ)
@btime repeval(hᵣ)
@btime repeval(h₂)
@btime repeval(h₂, A, y, λ)

#=
function l(t:: Type{<:Number}) where T
    x=ones(T, 10)

    return sum(map(sin, x))
end
=#