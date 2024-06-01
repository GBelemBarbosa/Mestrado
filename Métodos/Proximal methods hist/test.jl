using BenchmarkTools

function f(n:: Float64, s:: Vector{<:Number}, y:: Vector{<:Number})
    α=n/(s'y)
    α+=(α>10 || α<0.01)*(sqrt(n/(y'y))-α)
end 

function f₂(n:: Float64, s:: Vector{<:Number}, y:: Vector{<:Number})
    α=n/(s'y)
    if α>10 || α<0.01
        α=sqrt(n/(y'y))
    end
end 

function repeval(f:: Function)
    x=ones(100)
    for i in 1:10000
        res=f(10.0, x, 100 .*x)
    end
end

@btime repeval(f)
@btime repeval(f₂)

function g()
    x=1
    for i in 1:10000
        x=i
        x+=1
    end
end

function g₂()
    for i in 1:10000
        x=i
        x+=1
    end
end

@btime g()
@btime g₂()