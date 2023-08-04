using BenchmarkTools

f(x:: Vector{<:Number})=append!(x, [2*i for i=1:100])
function f₂(x:: Vector{<:Number})
    for i=1:100
        push!(x, 2*i)
    end
end

function repeval(f:: Function)
    for i in 1:10000
        res=f(ones(100))
    end
end

@btime repeval(f)
@btime repeval(f₂)