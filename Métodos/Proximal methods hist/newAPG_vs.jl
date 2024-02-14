using Plots
using LaTeXStrings

function newAPG_vs(F:: Function, g:: Function, ∂f:: Function, Tλ:: Function, Q:: Function, E:: Function, x₀:: Array{<:Number}, λ₁:: Number, μ₀:: Number, μ₁:: Number, c:: Number, δ:: Number, k_max:: Int64; ϵ=eps(), p=Inf) 
    x_, x=x₀, Tλ(λ₁, x₀, 0 .*x₀) #proxλ₁(x₀)
    Fx_, Fx=F(x_), F(x)
    ∂fx=∂f(x)
    ns=norm(x.-x_)^2
    λ=λ₁
    hist=[Fx_, Fx]
    
    k=1
    while true
        if ns<=c*(Fx_-Fx)
            y=x.+γ(k).*s
            ∂fy=∂f(y)
            x̂=Tλ(λ, y, ∂fy)
            Fx̂=F(x̂)
            
            if Fx̂<=Fx+min(Q(k), δ*(Fx_-Fx))
                Fy=F(y)
                x_, x=x, x̂
                xy=x.-y
                s=x.-x_
                nxy=xy'xy
                ns=s's
            else
                y=x
                Fy=Fx
                ∂fy=∂fx
                x_, x=x, Tλ(λ, y, ∂fy)
                xy=s=x.-x_
                nxy=ns=s's                
            end
        else
            y=x
            Fy=Fx
            ∂fy=∂fx
            x_, x=x, Tλ(λ, y, ∂fy)
            xy=s=x.-x_
            nxy=ns=s's
        end
        push!(hist, Fx)

        if norm(∂fy, p)<ϵ || k==k_max
            break
        end
        k+=1
        
        Fx_, Fx=Fx, F(x)
        ∂fx=∂f(x)
        aux=nxy/(λ*2*abs(g(x)Fy+g(y)-∂fx'xy-Fx))
        if 1>μ₀*aux
            λ=μ₁*aux
        else
            λ+=min(1, λₖ)*E(k)
        end
    end 

    return hist
end