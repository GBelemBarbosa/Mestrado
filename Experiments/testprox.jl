using BenchmarkTools

const λ=2

function proxhL(L:: Number, x:: Vector{<:Number}; λ=λ)
    λL=2*λ/L
    ωx=x.*x

    return [@inbounds x[i]*(ωx[i]>λL) for i=eachindex(x)]
end

function proxhL₂(L:: Number, x:: Vector{<:Number}; λ=λ)
    ωx=x.*x

    return [@inbounds x[i]*(ωx[i]>2*λ/L) for i=eachindex(x)]
end

function proxhL₃(L:: Number, x:: Vector{<:Number}; λ=λ)
    λL=2*λ/L

    return [@inbounds x[i]*(x[i]^2>λL) for i=eachindex(x)]
end

proxhL₄(L:: Number, x:: Vector{<:Number}; λ=λ)=[@inbounds x[i]*(x[i]^2>2*λ/L) for i=eachindex(x)]

proxhL₅(L:: Number, x:: Vector{<:Number}; λ=λ)=[@inbounds x[i]*(abs(x[i])>sqrt(2*λ/L)) for i=eachindex(x)]

proxhL₆(L:: Number, x:: Vector{<:Number}; λ=λ)=[@inbounds x[i]*(x[i]>sqrt(2*λ/L) || -x[i]>sqrt(2*λ/L)) for i=eachindex(x)]

function repeval(f:: Function)
    for i=1:1000
        res=f(5, rand(3000))
    end
end

#@btime repeval(proxhL)
#@btime repeval(proxhL₂)
#@btime repeval(proxhL₃)
@btime repeval(proxhL₄)
@btime repeval(proxhL₅)
@btime repeval(proxhL₆)