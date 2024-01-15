using Plots
using LaTeXStrings

function dual_projected_subgradient(f:: Function, g:: Function, oracle:: Function, γₖ:: Function, λ₀:: Array{<:Number}, k_max:: Int64; ϵ=eps(), p=Inf) 
    λ=λ₀
    x=oracle(λ)
    hist=[f(x)]
    
    k=1
    while true        
        gx=g(x)

        push!(hist, norm(max.(gx, 0), p))

        λ=max.(λ.+(γₖ(k, gx)/norm(gx, 2)).*gx, 0)
        x=oracle(λ)
        
        push!(hist, f(x))

        if norm(gx, p)<ϵ || k==k_max
            break
        end
        k+=1
    end 

    println(norm(max.(g(x), 0), p), " ", hist[end])
    x, scatter(eachindex(hist[begin:2:end]), [hist[begin:2:end], vcat(hist[2:2:end], 0)], 
                title=L"f(x^{(k)})"*" e "*L"||[g(x)]_+||_p",
                label=[L"f(x^{(k)})" L"||[g(x)]_+||_2"])
end