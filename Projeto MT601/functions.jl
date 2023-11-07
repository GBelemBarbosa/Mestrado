function quadratic(x:: Array{<:Number}; n=n)
    tot=0

    for i=1:n
        tot+=i*x[i]^2
    end

    return tot
end
∇quadratic(x:: Array{<:Number}; n=n)=[2*i*x[i] for i=1:n]
∇²quadratic(x:: Array{<:Number}; n=n)=diagm([2*i for i=1:n])
xₒₚₜquadratic(; n=n)=zeros(n)
fxₒₚₜquadratic(; n=n)=0

function rosenbrook(x:: Array{<:Number}; n=n)
    tot=0

    for i=1:Int64(n/2)
        tot+=10*(x[2*i]-x[2*i-1]^2)^2+(x[2*i-1]-1)^2
    end

    return tot
end
function ∇rosenbrook(x:: Array{<:Number}; n=n)
    ∇=Vector{Float64}(undef, n)

    for i=1:2:n
        ∇[i]=2*(x[i]-1)+40*x[i]*(x[i]^2-x[i+1])
    end
    for i=2:2:n
        ∇[i]=20*(x[i]-x[i-1]^2)
    end

    return ∇
end
function ∇²rosenbrook(x:: Array{<:Number}; n=n)
    ∇²=Array{Float64}(undef, n, n)

    for i=1:2:n
        ∇²[i, i]=2+40*(3*x[i]^2-x[i+1])
        aux=-40*x[i]
        ∇²[i, i+1]=aux
        ∇²[i+1, i]=aux
    end
    for i=2:2:n
        ∇²[i, i]=20
    end

    return ∇²
end
xₒₚₜrosenbrook(; n=n)=ones(n)
fxₒₚₜrosenbrook(; n=n)=0

function styblinsky_tang(x:: Array{<:Number}; n=n)
    tot=0

    for i=1:n
        tot+=x[i]^4-16*x[i]^2+5*x[i]
    end

    return tot
end
∇styblinsky_tang(x:: Array{<:Number}; n=n)=[4*x[i]^3-32*x[i]+5 for i=1:n]
∇²styblinsky_tang(x:: Array{<:Number}; n=n)=diagm([12*x[i]^2-32 for i=1:n])
xₒₚₜstyblinsky_tang(; n=n)=[-2.9035 for i=1:n]
fxₒₚₜstyblinsky_tang(; n=n)=-78.332*n

function rastrigin(x:: Array{<:Number}; n=n)
    tot=0

    for i=1:n
        tot+=x[i]^2-10*cos(2*pi*x[i])
    end

    return tot
end
∇rastrigin(x:: Array{<:Number}; n=n)=[2*x[i]+20*sin(2*pi*x[i])*pi for i=1:n]
∇²rastrigin(x:: Array{<:Number}; n=n)=diagm([2+40*cos(2*pi*x[i])*pi^2 for i=1:n])
xₒₚₜrastrigin(; n=n)=zeros(n)
fxₒₚₜrastrigin(; n=n)=-10*n

f_vec=[quadratic, rosenbrook, styblinsky_tang, rastrigin]
∇_vec=[∇quadratic, ∇rosenbrook, ∇styblinsky_tang, ∇rastrigin]
∇²_vec=[∇²quadratic, ∇²rosenbrook, ∇²styblinsky_tang, ∇²rastrigin]
xₒₚₜ_vec=[xₒₚₜquadratic, xₒₚₜrosenbrook, xₒₚₜstyblinsky_tang, xₒₚₜrastrigin]
fxₒₚₜ_vec=[fxₒₚₜquadratic, fxₒₚₜrosenbrook, fxₒₚₜstyblinsky_tang, fxₒₚₜrastrigin]