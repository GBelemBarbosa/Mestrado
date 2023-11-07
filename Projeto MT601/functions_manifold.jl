function quadratic(M, x:: Array{<:Number}; n=n)
    tot=0

    for i=1:n
        tot+=i*x[i]^2
    end

    return tot
end
∇quadratic(M, x:: Array{<:Number}; n=n)=[2*i*x[i] for i=1:n]

function rosenbrook(M, x:: Array{<:Number}; n=n)
    tot=0

    for i=1:Int64(n/2)
        tot+=10*(x[2*i]-x[2*i-1]^2)^2+(x[2*i-1]^2-1)^2
    end

    return tot
end
function ∇rosenbrook(M, x:: Array{<:Number}; n=n)
    ∇=Vector{Float64}(undef, n)

    for i=1:2:n
        ∇[i]=2*(x[i]-1)-40*x[i]*(x[i+1]-x[i]^2)
    end
    for i=2:2:n
        ∇[i]=20*(x[i]-x[i-1]^2)
    end

    return ∇
end

function styblinsky_tang(M, x:: Array{<:Number}; n=n)
    tot=0

    for i=1:n
        tot+=x[i]^4-16*x[i]^2+5*x[i]
    end

    return tot
end
∇styblinsky_tang(M, x:: Array{<:Number}; n=n)=[4*x[i]^3-32*x[i]+5 for i=1:n]

function rastrigin(M, x:: Array{<:Number}; n=n)
    tot=0

    for i=1:n
        tot+=x[i]^2-10*cos(2*pi*x[i])
    end

    return tot
end
∇rastrigin(x:: Array{<:Number}; n=n)=[2*x[i]+20*sin(2*pi*x[i])*pi for i=1:n]

f_vec=[quadratic, rosenbrook, styblinsky_tang, rastrigin]
∇_vec=[∇quadratic, ∇rosenbrook, ∇styblinsky_tang, ∇rastrigin]