using Plots
using LaTeXStrings

function dual_projected_subgradient(f:: Function, g:: Function, oracle:: Function, γₖ:: Function, λ₀:: Array{Number, N}, k_max:: Int64, ϵ:: Number) where {N}
    λ=λ₀
    x=oracle(λ)
    hist=[f(x)]
    
    for k=0:k_max        
        gx=g(x)

        push!(hist, norm(max.(gx, 0), 2))

        if norm(gx, Inf)<ϵ
            break
        end

        λ=max.(λ.+(γₖ(k, gx)/norm(gx, 2)).*gx, 0)
        x=oracle(λ)
        
        push!(hist, f(x))
    end 

    println(norm(max.(g(x), 0), 2), " ", f(x))
    scatter(eachindex(hist[begin:2:end]), [hist[begin:2:end], hist[2:2:end]], 
                title=L"f(x^{(k)})"*" e "*L"||[g(x)]_+||_2",
                label=[L"f(x^{(k)})" L"||[g(x)]_+||_2"])
end