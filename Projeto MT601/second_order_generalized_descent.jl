function SOGD(Hₖ:: Function, x₀:: Array{<:Number}, distxₒₚₜ:: Array{<:Number}, distfxₒₚₜ:: Array{<:Number}, tk_hist:: Array{<:Number}, converge:: Array{<:Number}, kϵ:: Array{Float64}; xₒₚₜ=xₒₚₜ, fxₒₚₜ=fxₒₚₜ, α=α, β=β, γ=γ, f=f, ∇f=∇f, tₖ=tₖ, resize=resize, fail=fail, t₀=t₀, H₀=H₀, k_max=k_max, s_max=s_max, ϵ_vec=ϵ_vec, p=p)
    d=x_=x=x₀
    t_=t=t₀
    H_=H=H₀
    fx_=fx=f(x)
    ∇fx=∇f(x)
    nd=n∇fx=norm(∇fx, p)
    threshold=0
    stop=length(ϵ_vec)

    for k=1:k_max
        lst_index=findall(n∇fx<i for i=ϵ_vec[threshold+1:end]).+threshold
        if !isempty(lst_index)
            threshold=lst_index[end]
            kϵ[lst_index].+=k/s_max
        end

        if threshold==stop #Stop criteria
            distx=norm(x.-xₒₚₜ, p)
            distxₒₚₜ[k:end].+=distx/s_max
            distfxₒₚₜ[k:end].+=(fx-fxₒₚₜ)/s_max
            if distx<=ϵ_vec[1]
                converge[1]+=1/s_max
            end

            return
        end

        while true 
            d=.-(H*∇fx)
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
            x=x_.+t.*d
            fx=f(x)

            if fx<=fx_+t*∇fxd #Armijo condition
                break
            end

            t=tₖ(k, t, x, d) #Step atualization
        end 

        distxₒₚₜ[k]+=norm(x.-xₒₚₜ, p)/s_max
        distfxₒₚₜ[k]+=(fx-fxₒₚₜ)/s_max
        tk_hist[k]+=t/s_max

        t_, t=t, resize(k, t_, t) #Step resize
        ∇fx_, ∇fx=∇fx, ∇f(x)
        n∇fx=norm(∇fx, p)
        H_, H=H, Hₖ(x_, x, ∇fx_, ∇fx, H_, H) #Hₖ atualization
        x_, fx_=x, fx
    end
    
    if norm(x.-xₒₚₜ, p)<=ϵ_vec[1]
        converge[1]+=1/s_max
    end
    kϵ[threshold+1:end].+=k_max/s_max
end