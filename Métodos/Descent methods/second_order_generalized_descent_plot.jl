using Plots
using LaTeXStrings

function FOGD(f:: Function, ∇f:: Function, tₖ:: Function, resize:: Function, Hₖ:: Function, fail:: Function, x₀:: Array{<:Number}, t₀:: Number, H₀:: Array{<:Number}, k_max:: Int64; α=0, β=Inf, γ=0, ϵ=eps, p=Inf)
    x_=x=x₀
    t_=t=t₀
    H_=H=H₀
    fx_=f(x)
    hist=[fx_]
    
    for k=0:k_max
        ∇fx=∇f(x)
        n∇fx=norm(∇fx, p)

        while true 
            d=.-H*∇f(x)
            nd=norm(d, p)
            ∇fxd=∇fx'd

            if ∇fxd<=-γ*n∇fx*nd #Descent direction condition 
                break 
            end

            H=fail(H_, H) #Procedure in case of failure
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
        
        if norm(x.-x_, p)<ϵ #Stop criteria
            break 
        end

        t_, t=t, resize(k, t_, t) #Step resize
        x_, fx_=x, fx
        ∇fx_=∇fx
        H_, H=H, Hₖ(x_, x, ∇fx_, ∇fx, H_, H) #Hₖ atualization
    end 

    println(norm(x.-x_, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end