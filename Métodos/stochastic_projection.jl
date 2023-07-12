using Random

function greedy_projection(PS:: Function, x₀:: Array{T, N}, m:: Int64, k_max:: Int64) where {T, N}
    x=x₀
    
    for k=0:k_max
        x=PS(x, rand(1:m)) #Projeta xₖ em Sᵢₖ, sendo iₖ randômico
    end 
end