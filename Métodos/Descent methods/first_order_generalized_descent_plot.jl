using Plots
using LaTeXStrings

function FOGD(f:: Function, ∇f:: Function, dₖ:: Function, tₖ:: Function, resize:: Function, x₀:: Array{<:Number}, t₀:: Number, k_max:: Int64; α=0, β=Inf, γ=0, ϵ=eps(), p=Inf)
    x_=x=x₀
    t_=t=t₀
    ∇fx_=∇fx=∇f(x)
    n∇fx=norm(∇fx, p)
    fx_=f(x)
    hist=[fx_]
    
    k=1
    while true
        while true 
            d=dₖ(k, x_, x, ∇fx_, ∇fx, d)
            nd=norm(d, p)
            ∇fxd=∇fx'd

            if ∇fxd<=-γ*n∇fx*nd #Descent direction condition 
                break 
            end
        end 

        n∇fx*=β/nd
        if 1<n∇fx #Direction size condition
            d.*=n∇fx
            ∇fxd*=n∇fx
        end
        
        ∇fxd*=α
        while true
            t=tₖ(k, t, x, d) #Step atualization
            x=x_.+t.*d
            fx=f(x)

            if fx<=fx_+t*∇fxd #Armijo condition
                break
            end
        end 

        push!(hist, fx)
        
        if norm(x.-x_, p) || k==k_max
            break
        end
        k+=1

        t_, t=t, resize(k, t_, t) #Step resize
        x_, fx_=x, fx
        ∇fx_, ∇fx=∇fx, ∇f(x)
        n∇fx=norm(∇fx, p)
    end 

    println(norm(x.-x_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end