using Plots
using LaTeXStrings

function FGCD(F:: Function, minimize:: Function, x₀:: Array{<:Number}; I₁=I₁, I₀=I₀, s=s)
    x=minimize(I₁(x₀))
    Fx=F(x)
    hist=Float64[]
    
    while true
        flag=true
        I₁x=I₁(x)
        g₀x=length(I₁x)
        options=[x]

        if !iszero(g₀x)
            append!(options, [minimize(filter(x->x!=i, I₁x)) for i=I₁x]) #xᵢ₋
        end
        if g₀x<=s
            I₀x=I₀(x)

            if g₀x<s
                append!(options, [minimize(vcat(I₁x, j)) for j=I₀x]) #xⱼ₊
            end
            if !iszero(g₀x)
                for i=I₁x
                    append!(options, [minimize(vcat(filter(x->x!=i, I₁x), j)) for j=I₀x]) #xᵢⱼ
                end
            end
        end
        for xᵢ=options
            Fxᵢ=F(xᵢ)
            if Fx>Fxᵢ
                flag=false
                x=xᵢ
                Fx=Fxᵢ
            end
        end
        
        if flag
            break
        end

        push!(hist, Fx)
    end 

    println(Fx)
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end