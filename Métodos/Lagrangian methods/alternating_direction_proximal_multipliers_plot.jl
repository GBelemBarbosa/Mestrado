using Plots
using LaTeXStrings

function AD_PMM(H:: Function, A:: Array{<:Number, M}, B:: Array{<:Number, M}, c:: Array{<:Number, N}, minimize_x:: Function, minimize_z:: Function, ρ:: Number, y₀:: Array{<:Number, N}, k_max:: Int64; ϵ=eps, p=Inf) where {M, N}
    y=y₀
    x=A'*y
    z=B'*y
    Bz=B*z
    hist=[hist[end]]
    
    for k=0:k_max        
        aux=y./ρ.-c
        x_, x=x, minimize_x(x, Bz.+aux)
        Ax=A*x
        z_, z=z, minimize_z(z, Ax.+aux)

        push!(hist, H(x))

        if max(norm(x.-x_, p), norm(z.-z_, p))<ϵ
            break
        end

        y.+=ρ.*(Ax.+By.-c)
    end 

    println(hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"H(x^{(k)})",
                label=false)
end