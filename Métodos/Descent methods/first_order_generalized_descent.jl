using Plots
using LaTeXStrings

function FOGD(f:: Function, ∇f:: Function, dₖ:: Function, step:: Function, x₀:: Array{<:Number}, t₀:: Number, k_max:: Int64; α=0, β=Inf, γ=0, ϵ=eps, p=Inf)
    x=x₀
    t=t₀
    hist=[f(x)]
    
    for k=0:k_max
        ∇fx=∇f(x)
        n∇fx=norm(∇fx, p)

        while true 
            d=dₖ(x, d, ∇fx)
            nd=norm(d, p)

            if ∇fx'd<=γ*n∇fx*nd #Descent direction condition 
                break 
            end
        end 

        if nd<β*n∇fx #Direction size condition
            d*=β*n∇fx/nd
            nd=β*n∇fx
        end
        
        x, t=step(k, x, t, d, α) #Backtracking or another armijo-guarantee procedure

        push!(hist, f(x))

        if nd<ϵ #Stop criteria
            break 
        end
    end 

    println(norm(d, p), " ", hist[end])
    x, scatter(eachindex(hist), hist, 
                title=L"f(x^{(k)})",
                label=false)
end