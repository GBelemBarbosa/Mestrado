ω(x:: Vector{<:Number})=x.*x

T(ωx:: Vector{<:Number})=sortperm(ωx, rev=true)

function proxl0L(L:: Number, x:: Vector{<:Number}, λ:: Number; ω=ω)
    λL=2*λ/L
    ωx=ω(x)

    return [x[i]*(ωx[i]>λL) for i=eachindex(x)]
end