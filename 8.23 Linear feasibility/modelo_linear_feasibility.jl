include("variables_linear_feasibility.jl")

PS₁(x:: Vector{<:Number})=max.(x, 0) #Projeção em R₊ⁿ

d₁(x:: Vector{<:Number})=abs(minimum(x)) #Distância à R₊ⁿ

PS₂(x:: Vector{<:Number}; A=A, b=b)=x.-A'*((A*A')\(A*x.-b)) #Projeção em {x: Ax=b}

include("../Métodos/Projected methods/alternating_projection_plot.jl")

p₁=alternating_projection(PS₁, PS₂, d₁, x₀, k_max, ϵ)

PS₁(x:: Vector{<:Number}; A=A, b=b)=x.-A'*((A*A')\(A*x.-b)) #Projeção em {x: Ax=b}

d₁(x:: Vector{<:Number})=norm(abs.(A*x-b), Inf) #Distância à {x: Ax=b} em termos de A*x-b

PS₂(x:: Vector{<:Number})=max.(x, 0) #Projeção em R₊ⁿ

p₂=alternating_projection(PS₁, PS₂, d₁, x₀, k_max, ϵ)

plot(p₁, p₂)