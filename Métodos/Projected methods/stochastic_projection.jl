using Random

function sthocastic_projection(PS:: Function, x₀:: Array{Number, N}, m:: Int64, k_max:: Int64) where {N}
    x=x₀
    
    for k=0:k_max
        x=PS(x, rand(1:m)) #Projeta xₖ em Sᵢₖ, sendo iₖ randômico (pode ser implementado algo que garante dois índices consecutivos diferentes)
    end 
end