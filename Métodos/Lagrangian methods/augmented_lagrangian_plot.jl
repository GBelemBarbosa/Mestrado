using Plots
using LaTeXStrings

function augmented_lagrangian(H:: Function, A:: Array{<:Number, M}, B:: Array{<:Number, M}, c:: Array{<:Number, N}, minimize_Lρ:: Function, ρ:: Number, y₀:: Array{<:Number, N}, k_max:: Int64; ϵ=eps, p=Inf) where {M, N}
    y=y₀
    x=A'*y
    z=B'*y
    hist=[H(x)]
    
    for k=0:k_max        
        x_, z_, x, z=x, z, minimize_Lρ(y)

        push!(hist, H(x))

        if max(norm(x.-x_, p), norm(z.-z_, p))<ϵ
            break
        end

        y.+=ρ.*(A*x.+B*y.-c)
    end 

    println(hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"H(x^{(k)})",
                label=false)
end