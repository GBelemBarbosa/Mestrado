g(x:: Vector{Number}, m:: Int64, group_indexs:: Vector{UnitRange{Int64}})=[!iszero(x[group_indexs[i]]) for i=1:m]

g₀(x:: Vector{Number}; g=g)=sum(g(x))

A(T:: Int64, group_indexs:: Vector{UnitRange{Int64}})=group_indexs[T]

A(T:: Vector{Int64}, group_indexs:: Vector{UnitRange{Int64}})=reduce(vcat, group_indexs[i] for i=T)

δCₛB(x:: Vector{Number}, s:: Int64; B=B)=0+(g₀(x)>s+!all(B(x)))*Inf

h(x:: Vector{Number}, λ:: Number; δCₛB=δCₛB, g₀=g₀)=λ*g₀(x)+δCₛB(x)

function UATPBTAT(x:: Vector{Number}, j:: Int64, n:: Int64, group_indexs:: Vector{UnitRange{Int64}}; A=A)
    y=zeros(Float64, n)
    y[A(j, group_indexs)]=PDⱼ(x[A(j, group_indexs)], j)

    return y
end

function UATPBTAT(x:: Vector{Number}, T:: Vector{Int64}, n:: Int64, group_indexs:: Vector{UnitRange{Int64}}; A=A)
    y=zeros(Float64, n)
    for j=T
        y[A(j, group_indexs)]=PDⱼ(x[A(j, group_indexs)], j)
    end

    return y
end

PBTAT(x:: Vector{Number}, j:: Int64, group_indexs:: Vector{UnitRange{Int64}}; A=A)=PDⱼ(x[A(j, group_indexs)], j)

PBTAT(x:: Vector{Number}, T:: Vector{Int64}, group_indexs:: Vector{UnitRange{Int64}}; A=A)=reduce(vcat, PDⱼ(x[A(j, group_indexs)], j) for j=T)

ω(x:: Vector{Number}, m:: Int64, group_indexs:: Vector{UnitRange{Int64}}; A=A)=[norm(x[A(j, group_indexs)], 2)^2-dDⱼ(x[A(j, group_indexs)], j)^2 for j=1:m]

ωₛ(x:: Vector{Number}, s=Inf64; ω=ω)=partialsort(ω(x), s, rev=true)

function Sₛ(x:: Vector{Number}, s=Inf64; ω=ω)
    ωx=ω(x)
    uωx=unique(ωx)
    uωxₛ=partialsort(uωx, min(s, length(uωx)), rev=true)

    return findall(x -> x>=uωxₛ, ωx)
end

I₁(x:: Vector{Number}, group_indexs:: Vector{UnitRange{Int64}}, m:: Int64)=findall(!iszero, g(x, m, group_indexs))

I₀(x:: Vector{Number}, group_indexs:: Vector{UnitRange{Int64}}, m:: Int64)=findall(iszero, g(x, m, group_indexs))

function I₊(x:: Vector{Number}, group_indexs:: Vector{UnitRange{Int64}}, m:: Int64, s:: Int64, λ:: Int64)
    ωx=ω(x, m:: Int64, group_indexs)
    ωxₛ=partialsort(ωx, s, rev=true)

    return findall(x->x>max(ωxₛ, 2*λ), ωx)
end

function Iq(x:: Vector{Number}, s:: Int64, λ:: Int64)
    ωx=ω(x, m:: Int64, group_indexs)
    ωxₛ=partialsort(ωx, s, rev=true)

    return findall(isequal(max(ωxₛ, 2*λ)), ωx)
end

T(ωx:: Vector{Number}, s:: Int64)=partialsortperm(ωx, 1:s, rev=true)

function proxhL(L:: Number, x:: Vector{Number}, s:: Int64, λ:: Int64, n:: Int64, group_indexs:: Vector{UnitRange{Int64}}; T=T, ω=ω)
    ωx=ω(x)
    Tωx=T(ωx)

    return UATPBTAT(x, Tωx[1:searchsortedlast(ωx[Tωx], 2*λ/L, rev=true, lt=<=)], n, group_indexs)
end

TL(L:: Number, x:: Vector{Number}; ∇f=∇f)=x.-∇f(x)./L