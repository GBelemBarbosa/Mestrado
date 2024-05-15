function IBPG(∂f:: Function, ℘hλg:: Function, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps(), p=Inf) 
    z=x_=x=x₀
    αₖ=βₖ=λₖ=0
    
    k=1
    while true
        s=(x.-x_)
        x_=copy(x)
        z_, z=z, x.+βₖ.*s
        ∂fz_, ∂fz=∂fz, ∂f(z)
        αₖ, βₖ, λₖ, x=℘hλg(αₖ, βₖ, λₖ, k, s, x, z_, z, ∂fz_, ∂fz) #Atualização

        if norm(∂fz, p)<ϵ || k==k_max
            break
        end
        k+=1
    end 

    return x
end