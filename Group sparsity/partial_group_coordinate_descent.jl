using Plots
using LaTeXStrings

function PGCD(F:: Function, minimize:: Function, L:: Number, x₀:: Array{<:Number}, k_max:: Int64; I₁=I₁, ω=ω, g=g, iₓₗ=iₓₗ, jₓₗ=jₓₗ)
    x=minimize(I₁(x₀))
    Fx=F(x)
    hist=Float64[]
    
    for k=0:k_max
        flag=true
        I₁x=I₁(x)
        g₀x=length(I₁x)
        ωTLx=ω(TL(L, x))
        gx=g(x)
        i=iₓₗ(ωTLx, gx)
        options=[x]

        if !iszero(g₀x)
            push!(options, minimize(filter(x->x!=i, I₁x))) #xᵢ₋
        end
        if g₀x<=s
            j=jₓₗ(ωTLx, gx)

            if g₀x<s
                push!(options, minimize(vcat(I₁x, j))) #xⱼ₊
            end
            if !iszero(g₀x)
                push!(options, minimize(vcat(filter(x->x!=i, I₁x), j))) #xᵢⱼ
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

    println(F(x))
    x, scatter(eachindex(hist), hist, 
                title=L"F(x^{(k)})",
                label=false)
end