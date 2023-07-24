n=[3, 8, 15]
m=length(n)

pushfirst!(n, 0)
group_indexs=[n[i]+1:n[i+1] for i=1:m]

g(x:: Vector{<:Number}; group_indexs=group_indexs, m=m)=[!iszero(x[group_indexs[i]]) for i=1:m]

g₀(x:: Vector{<:Number}; g=g)=sum(g(x))

A(T:: Int64; group_indexs=group_indexs)=group_indexs[T]

A(T:: Vector{Int64}; group_indexs=group_indexs)=[group_indexs[i] for i=T]

s=2

B(x:: Vector{<:Number})=(all(x[A(1)].<=0), all(x[A(2)].>=0), x[A(3)][3]>=-2)

δCₛB(x:: Vector{<:Number}; B=B, s=s)=0+(g₀(x)>s+!all(B(x)))*Inf