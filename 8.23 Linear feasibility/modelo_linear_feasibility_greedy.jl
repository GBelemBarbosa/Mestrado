include("variables_linear_feasibility.jl")

function PS(x:: Vector{Float64}, i:: Int64; A=A, b=b, m=m) #Projeções
    if i>m
        return max.(x, 0)
    else
        return x.+((b[i]-A[i, :]'*x)/norm(A[i, :])^2).*A[i, :]
    end
end

function d(x:: Vector{Float64}; A=A, b=b, m=m) #Distância máxima e índice da mesma
    d_max=0

    for i=1:m
        dᵢ=abs(A[i, :]'*x-b[i])/norm(A[i, :])
        if dᵢ>d_max
            d_max=dᵢ
            i_max=i
        end
    end
    dᵢ=norm(x.-max.(x, 0))
    if dᵢ>d_max
        d_max=dᵢ
        i_max=m+1
    end

    return d_max, i_max
end

include("../Métodos/Projected methods/greedy_projection_plot.jl")

p=greedy_projection(PS, d, x₀, k_max, ϵ)