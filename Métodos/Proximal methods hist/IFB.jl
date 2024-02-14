using Plots
using LaTeXStrings

function IFB(∂g:: Function, argmin:: Function, x₀:: Array{<:Number}, n_max:: Int64; ϵ=eps(), p=Inf) 
    x_=x=x₀
    αₙ=βₙ=0
    
    n=1
    while true
        ∂gx_, ∂gx=∂gx, ∂g(x)
        x_, αₙ, βₙ, x, x_=x, argmin(αₙ, βₙ, n, x_, x, ∂gx_, ∂gx) #Backtracking mais atualização

        if norm(∂gx, p)<ϵ || n==n_max
            break
        end
        n+=1
    end 

    return x
end