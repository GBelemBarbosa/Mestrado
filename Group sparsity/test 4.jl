using BenchmarkTools

f(x:: Vector{<:Number}, m:: Int64)=[x[i] for i=1:m]

f₂(x:: Vector{<:Number}, m:: Int64)=[@inbounds x[i] for i=1:m]

f(x:: Vector{<:Number})=[x[i] for i=eachindex(x)]

indexs=1:100

h(x:: Vector{<:Number}, m:: Int64; indexs=indexs)=[x[indexs[i]] for i=1:m]

h₂(x:: Vector{<:Number}, m:: Int64; indexs=indexs)=[@inbounds x[indexs[i]] for i=1:m]

h(x:: Vector{<:Number}; indexs=indexs)=[x[indexs[i]] for i=eachindex(indexs)]

h₃(x:: Vector{<:Number}; indexs=indexs)=[@inbounds x[indexs[i]] for i=eachindex(indexs)]

function repeval(f:: Function)
    x=ones(100)
    for i in 1:10000
        res=f(x)
    end
end

function repeval(f:: Function, m:: Int64)
    x=ones(100)
    for i in 1:10000
        res=f(x, m)
    end
end

@btime repeval(f)
@btime repeval(f₂, 100)
@btime repeval(f, 100)
@btime repeval(h)
@btime repeval(h₃)
@btime repeval(h₂, 100)
@btime repeval(h, 100)