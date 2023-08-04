using BenchmarkTools
using LinearAlgebra

f(x:: Vector{<:Number})=dot(ones(100), x)
f₂(x:: Vector{<:Number})=x'*ones(100)
f₃(x:: Vector{<:Number})=x'ones(100)

function repeval(f:: Function)
    for i in 1:10000
        res=f(ones(100))
    end
end

@btime repeval(f)
@btime repeval(f₂)
@btime repeval(f₃)