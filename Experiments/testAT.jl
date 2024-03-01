using BenchmarkTools

A=ones(5000, 5000)

function repeval(; A=A)
    f(x:: Vector{<:Number}; A=A)=A'x

    for i=1:1000
        res=f(ones(5000))
    end
end

function repeval₂(; A=A)
    f(x:: Vector{<:Number}; A=A')=A*x
    
    for i=1:1000
        res=f(ones(5000))
    end
end

function repeval₃(; A=A)
    AT=A'
    f(x:: Vector{<:Number}; A=AT)=A*x

    for i=1:1000
        res=f(ones(5000))
    end
end

@btime repeval()
@btime repeval₂()
@btime repeval₃()

f(x:: Vector{<:Number}; A=A)=A'x

f₂(x:: Vector{<:Number}; A=A)=A'*x

f₃(x:: Vector{<:Number}; A=A)=(x'*A)'

function repeval(f:: Function; A=A)
    for i=1:1000
        res=f(ones(5000))
    end
end

@btime repeval(f)
@btime repeval(f₂)
@btime repeval(f₃)