using BenchmarkTools

f(x:: Vector{<:Number})=[x[i]^2 for i=eachindex(x)]

f₂(x:: Vector{<:Number})=x.*x

function repeval(f:: Function)
    x=ones(100)
    for i in 1:10000
        res=f(x)
    end
end

@btime repeval(f)
@btime repeval(f₂)