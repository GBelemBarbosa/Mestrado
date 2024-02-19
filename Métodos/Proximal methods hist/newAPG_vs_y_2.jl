using Plots
using LaTeXStrings

function newAPG_vs_y_2(F:: Function, g:: Function, ∂f:: Function, Tλ:: Function, γ:: Function, Q:: Function, E:: Function, x₀:: Array{<:Number}, λ₁:: Number, μ₀:: Number, μ₁:: Number, c:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    x_, x=x₀, Tλ(λ₁, x₀, 0 .*x₀) # proxλ₁(x₀)
    Fx_, Fx=F(x_), F(x)
    s=x.-x_
    ns=norm(s)^2
    λ=λ₁
    t=1.0
    hist=[Fx_, Fx]
    histnψ=[]
    ls=0
    
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
                Fx_, Fx=Fx, F(x)
                xy=x.-y
                nxy=xy'xy
                s=x.-x_
                ns=s's
            else
                ỹ=x
                Fy=Fx
                ∂fỹ=∂f(y)
                x_, x=x, Tλ(λ, ỹ, ∂fỹ)

                ls+=1

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
                    y=ỹ
                    ∂fy=∂fỹ
                    xy=s=x.-x_
                    nxy=ns=s's 
                end             
            end
        else
            y=x
            Fy=Fx
            ∂fy=∂f(y)
            x_, x=x, Tλ(λ, y, ∂fy)

            ls+=1

            Fx_, Fx=Fx, F(x)
            xy=s=x.-x_
            nxy=ns=s's
        end

        push!(hist, Fx)
        nψ=norm(∂f(x).-∂fy.+(y.-x)./λ, p)
        push!(nψ, histnψ)

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

    return x, hist, histnψ, ls
end