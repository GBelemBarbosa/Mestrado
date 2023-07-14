using BenchmarkTools

f(x:: Vector{Float64})=x-x/2
f₂(x:: Vector{Float64})=x.-x./2

g(x:: Vector{Float64})=1 .-x./2
g₂(x:: Vector{Float64})=ones(100).-x./2

function repeval(f:: Function)
    for i in 1:10000
        res = f(ones(100))
    end
end

@btime repeval(f)
@btime repeval(f₂)

@btime repeval(g)
@btime repeval(g₂)

h(x:: Number)=1-x/2
h₂(x:: T) where T=1-x/2

function repeval(f:: Function)
    for i in 1:10000
        res = f(1)
    end
end

@btime repeval(h)
@btime repeval(h₂)