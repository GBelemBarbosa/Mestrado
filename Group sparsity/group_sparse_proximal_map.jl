include("problem_formulation.jl")

g(x:: Vector{<:Number}; group_indexs=group_indexs, m=m)=[!iszero(x[group_indexs[i]]) for i=1:m]

g₀(x:: Vector{<:Number}; g=g)=sum(g(x))

A(T:: Int64; group_indexs=group_indexs)=group_indexs[T]

A(T:: Vector{Int64}; group_indexs=group_indexs)=reduce(vcat, group_indexs[i] for i=T)

δCₛB(x:: Vector{<:Number}; B=B, s=s)=0+(g₀(x)>s+!all(B(x)))*Inf

function UATPBTAT(x:: Vector{<:Number}, j:: Int64; A=A, m=m, n=n)
    y=zeros(Float64, n[m+1])
    y[A(j)]=PDⱼ(x[A(j)], j)

    return y
end

function UATPBTAT(x:: Vector{<:Number}, T:: Vector{Int64}; A=A, m=m, n=n)
    y=zeros(Float64, n[m+1])
    for j=T
        y[A(j)]=PDⱼ(x[A(j)], j)
    end

    return y
end

PBTAT(x:: Vector{<:Number}, j:: Int64; A=A, m=m, n=n)=PDⱼ(x[A(j)], j)

PBTAT(x:: Vector{<:Number}, T:: Vector{Int64}; A=A, m=m, n=n)=reduce(vcat, PDⱼ(x[A(j)], j) for j=T)

ω(x:: Vector{<:Number}; A=A, m=m)=[norm(x[A(j)], 2)^2-dDⱼ(x[A(j)], j)^2 for j=1:m]

ωₛ(x:: Vector{<:Number}; ω=ω, s=s)=partialsort(ω(x), s, rev=true)

function Sₛ(x:: Vector{<:Number}; ω=ω, s=s)
    ωx=ω(x)
    uωx=unique(ωx)
    uωxₛ=partialsort(uωx, min(s, length(uωx)), rev=true)

    return findall(x->x>=uωxₛ, ωx)
end

I1(x:: Vector{<:Number})=findall(!iszero, g(x))

function I₊(x:: Vector{<:Number}; s=s, λ=λ)
    ωx=ω(x)
    ωxₛ=partialsort(ωx, s, rev=true)

    return findall(x->x>max(ωxₛ, 2*λ), ωx)
end

function Iq(x:: Vector{<:Number}; s=s, λ=λ)
    ωx=ω(x)
    ωxₛ=partialsort(ωx, s, rev=true)

    return findall(isequal(max(ωxₛ, 2*λ)), ωx)
end

T(ωx:: Vector{<:Number}; s=s)=partialsortperm(ωx, 1:s, rev=true)

function proxₕ(x:: Vector{<:Number}; ω=ω, s=s, λ=λ)
    ωx=ω(x)
    R=[i for i=T(ωx) if ωx[i]>2*λ]

    return UATPBTAT(x, R[1:min(s, length(R))])
end

TL(L:: Number, x:: Vector{<:Number}; ∇f=∇f)=x.-∇f(x)./L

iₓₗ(ωTLx:: Vector{<:Number}, gx:: Vector{Bool}; TL=TL)=argmin(ωTLx.*gx)

jₓₗ(ωTLx:: Vector{<:Number}, gx:: Vector{Bool}; TL=TL)=argmax(ωTLx.*.!gx)

function ijₓₗ(L:: Number, x:: Vector{<:Number}; ω=ω, TL=TL, g=g)
    ωTLx=TL(L, ω(x))
    gx=g(x)

    return iₓₗ(ωTLx, gx), iₓₗ(ωTLx, gx)
end