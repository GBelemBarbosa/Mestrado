using BenchmarkTools

function f(x:: Vector{Float64})
    x_=ones(100)
    x_, x[1:20]=x, zeros(20)
end

function f₂(x:: Vector{Float64})
    x_=ones(100)
    x_[1:20], x[1:20]=x[1:20], zeros(20)
end

function repeval(f:: Function)
    for i in 1:10000
        res = f(ones(100))
    end
end

@btime repeval(f)
@btime repeval(f₂)