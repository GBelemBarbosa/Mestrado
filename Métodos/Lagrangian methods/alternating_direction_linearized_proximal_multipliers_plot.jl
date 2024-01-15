using Plots
using LaTeXStrings

function AD_LPMM(H:: Function, A:: Array{<:Number, M}, B:: Array{<:Number, M}, c:: Array{<:Number, N}, proxh₁α:: Function, proxh₂β:: Function, ρ:: Number, α:: Number, β:: Number, y₀:: Array{<:Number, N}, k_max:: Int64; ϵ=eps(), p=Inf) where {M, N}
    y=y₀
    x=A'*y
    Ax=A*x
    z=B'*y
    Bz=B*z
    Aux=Ax.+Bz
    hist=[H(x)]
    
    k=1
    while true 
        aux=y./ρ.-c       
        x_, x=x, proxh₁α(x.-(ρ/α).*A'*(Aux.+aux))
        Ax=A*x
        z_, z=z, proxh₂β(z.-(ρ/β).*B'*(Ax.+Bz.+aux))
        Bz=B*z

        push!(hist, H(x))

        if max(norm(x.-x_, p), norm(z.-z_, p))<ϵ || k==k_max
            break
        end
        k+=1

        Aux=Ax.+Bz
        y.+=ρ.*(Aux.-c)
    end 

    println(hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"H(x^{(k)})",
                label=false)
end