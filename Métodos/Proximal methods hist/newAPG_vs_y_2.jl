using Plots
using LaTeXStrings

function newAPG_vs_y_2(F:: Function, g:: Function, ∂f:: Function, Tλ:: Function, γ:: Function, Q:: Function, E:: Function, x₀:: Array{<:Number}, λ₁:: Number, μ₀:: Number, μ₁:: Number, c:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    t1=time()
    ∂fx=∂f(x₀)
    start=time()
    maybe=t1-start
    λ=λ₁
    t=1.0
    histnψ=Tuple{Float64, Float64}[]
    ls=0
    Fx_=F(x₀)
    histF=[(time()-start, Fx_)]
    x_, x=x₀, Tλ(λ₁, x₀, 0 .*x₀) # proxλ₁(x₀)
    Fx=F(x)
    s=x.-x_
    ns=norm(s)^2
    push!(histF, (time()-start, Fx))
    
    k=1
    k_max-=1
    while true
        if ns<=c*(Fx_-Fx)
            t, γₖ=γ(k, t)
            y=x.+γₖ.*s
            ∂fy=∂f(y)
            x̂=Tλ(λ, y, ∂fy)

            ls+=1

            Fx̂=F(x̂)
            
            if Fx̂<=Fx+min(Q(k), δ*(Fx_-Fx))
                Fy=F(y)
                x_, x=x, x̂
                Fx_, Fx=Fx, Fx̂
                xy=x.-y
                nxy=xy'xy
                s=x.-x_
                ns=s's
            else
                x_, x=x, Tλ(λ, x, ∂fx)

                ls+=1
                start+=maybe

                Fx_, Fx=Fx, F(x)
                if Fx>Fx̂
                    Fy=F(y)
                    x=x̂
                    Fx=Fx̂
                    xy=x.-y
                    nxy=xy'xy
                    s=x.-x_
                    ns=s's
                else
                    Fy=Fx_
                    y=x
                    ∂fy=∂fx
                    xy=s=x.-x_
                    nxy=ns=s's 
                end             
            end
        else
            y=x
            Fy=Fx
            ∂fy=∂fx
            x_, x=x, Tλ(λ, y, ∂fy)

            ls+=1
            start+=maybe

            Fx_, Fx=Fx, F(x)
            xy=s=x.-x_
            nxy=ns=s's
        end
        t1=time()
        ∂fx=∂f(x)
        maybe=t1-time()

        t1=time()
        elapsed=t1-start
        nψ=norm(∂fx.-∂fy.+(y.-x)./λ, p)
        push!(histF, (elapsed, Fx))
        push!(histnψ, (elapsed, nψ))
        start+=time()-t1

        if nψ<ϵ || k==k_max
            break
        end
        k+=1
        
        aux=nxy/(2*abs(Fy-g(y)-Fx+g(x)+∂fy'xy))
        if 1>μ₀*aux/λ
            λ=μ₁*aux
        else
            λ+=min(1, λ)*E(k)
        end
    end 

    return x, histF, histnψ, ls
end