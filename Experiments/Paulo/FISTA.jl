using Plots
using LaTeXStrings

function FISTA(F:: Function, ∂f:: Function, proxhL:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    Fx=F(x₀) # Shouldn't count against FISTA time (because isn't necessary to the method)
    start=time()
    y=x_=x=x₀
    t=1
    L=s
    histnψ=Tuple{Float64, Float64}[]
    histF=[(time()-start, Fx)]
    
    k=1
    while true
        ∂fy=∂f(y)
        x=proxhL(L, y.-∂fy./L) 

        t1=time()
        elapsed=t1-start
        nψ=norm(∂f(x).-∂fy.+(y.-x).*L, p)
        push!(histF, (elapsed, F(x)))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1

        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_)
        x_=x
    end 

    return x, histF, histnψ
end