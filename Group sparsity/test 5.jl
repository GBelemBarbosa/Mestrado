using BenchmarkTools

function f(x:: Vector{<:Number}, λL:: Number)
    for i=eachindex(x)
        if x[i]<λL
            if x[i]<-λL
                x[i]+=λL
            else 
                x[i]=0.0
            end
        else
            x[i]-=λL
        end
    end

    return x
end

function f₂(x:: Vector{<:Number}, λL:: Number)
    for i=eachindex(x)
        if x[i]>λL
            x[i]-=λL
        elseif x[i]<-λL
            x[i]+=λL
        else
            x[i]=0.0
        end
    end

    return x
end

function f₃(x:: Vector{<:Number}, λL:: Number)
    aux=zeros(Float64, length(x))

    for i=eachindex(x)
        if x[i]<λL
            if x[i]<-λL
                aux[i]=x[i]+λL
            end
        else
            aux[i]=x[i]-λL
        end
    end

    return x
end

f₄(x:: Vector{<:Number}, λL:: Number)=[0.0+(x[i]<-λL)*(x[i]+λL)+(x[i]>λL)*(x[i]-λL) for i=eachindex(x)]

f₅(x:: Vector{<:Number}, λL:: Number)=[x[i]+(x[i]<-λL)*λL-(x[i]>λL)*λL-(abs(x[i])<λL)*x[i] for i=eachindex(x)]

function repeval(f:: Function)
    for i in 1:10000
        res=f(randn(100), 0.1)
    end
end

@btime repeval(f)
@btime repeval(f₂)
@btime repeval(f₃)
@btime repeval(f₄)
@btime repeval(f₅)
