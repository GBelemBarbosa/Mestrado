using BenchmarkTools

A=ones(500, 500)

function repeval(; A=A)
    f(x:: Vector{<:Number}; A=A)=A'x

    for i in 1:100
        res=f(ones(500))
    end
end

function repeval₂(; A=A)
    f(x:: Vector{<:Number}; A=A')=A*x
    
    for i in 1:100
        res=f(ones(500))
    end
end

function repeval₃(; A=A)
    AT=A'
    f(x:: Vector{<:Number}; A=AT)=A*x

    for i in 1:100
        res=f(ones(500))
    end
end

@btime repeval()
@btime repeval₂()
@btime repeval₃()