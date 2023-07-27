using LinearAlgebra

n=[3, 8, 15]
m=length(n)
λ=1
s=2

pushfirst!(n, 0)
group_indexs=[n[i]+1:n[i+1] for i=1:m]

B(x:: Vector{<:Number})=(all(x[A(1)].<=0), all(x[A(2)].>=0), x[A(3)][3]>=-2)

x=randn(Float64, n[m+1])

function PDⱼ(xGⱼ:: Vector{<:Number}, j:: Int64; A=A)
    if j==1
        return min.(xGⱼ, 0)
    elseif j==2
        return max.(xGⱼ, 0)
    end
    xGⱼ[3]=max(xGⱼ[3], -2)

    return xGⱼ
end

function PDⱼAj(x:: Vector{<:Number}, j:: Int64; A=A)
    if j==1
        return min.(x[A(1)], 0)
    elseif j==2
        return max.(x[A(2)], 0)
    end
    xG₃=x[A(3)]
    xG₃[3]=max(xG₃[3], -2)

    return xG₃
end

function PB(x:: Vector{<:Number}; A=A)
    xG3=x[A(3)]
    xG3[3]=max(xG3[3], -2)

    return vcat(min.(x[A(1)], 0), max.(x[A(2)], 0), xG3)
end

function dDⱼ(xGⱼ:: Vector{<:Number}, j:: Int64; A=A)
    if j==1
        return norm(max.(xGⱼ, 0), 2)
    elseif j==2
        return norm(min.(xGⱼ, 0), 2)
    end

    return max(-xGⱼ[3], 2)-2
end

function dDⱼAj(x:: Vector{<:Number}, j:: Int64; A=A)
    if j==1
        return norm(max.(x[A(1)], 0), 2)
    elseif j==2
        return norm(min.(x[A(2)], 0), 2)
    end

    return max(-x[A(3)][3], 2)-2
end

dB(x:: Vector{<:Number}; A=A)=[norm(max.(x[A(1)], 0), 2), norm(min.(x[A(2)], 0), 2), max(-x[A(3)][3], 2)-2]

function ωⱼ(x:: Vector{<:Number}, j:: Int64; A=A, dDⱼ=dDⱼ)
    xGⱼ=x[A(j)]
    if j==1
        return norm(xGⱼ, 2)^2-dDⱼ(xGⱼ, 1)^2
    elseif j==2
        return norm(xGⱼ, 2)^2-dDⱼ(xGⱼ, 2)^2
    end

    return norm(xGⱼ, 2)^2-dDⱼ(xGⱼ, 3)^2
end

∇f(x:: Vector{<:Number})=x.-1