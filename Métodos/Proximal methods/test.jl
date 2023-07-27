using BenchmarkTools

function f(x:: Vector{<:Number})
    x[1:20].+=1
    b=x[1:20]
end

function f₂(x:: Vector{<:Number})
    aux=x[1:20].+=1
    b=aux
end

function repeval(f:: Function)
    for i in 1:10000
        res = f(ones(100))
    end
end

@btime repeval(f)
@btime repeval(f₂)

function g(x:: Vector{<:Number})
    b=ones(20)
    x[1:20].+=1
    b.+=x[1:20]
end

function g₂(x:: Vector{<:Number})
    b=ones(20)
    aux=x[1:20].+=1
    b.+=aux
end

@btime repeval(g)
@btime repeval(g₂)

function h(x:: Vector{<:Number})
    a, b=x, x
end

function h₂(x:: Vector{<:Number})
    a=b=x
end

@btime repeval(h)
@btime repeval(h₂)