include("variables_Huber.jl")

hη(w:: Float64; η=η)=(abs(w)>η)*η*(abs(w)-η/2)+!(abs(w)>η)*w^2/2
f(x:: Vector{Float64}; A=A, y=y, λ=λ, hη=hη)=sum(hη.(y.-A*x))+λ*norm(x)^2/2

∂hη(w:: Float64; η=η)=(abs(w)>η)*η*sign(w)+!(abs(w)>η)*w
∇f(x:: Vector{Float64}; A=A, y=y, λ=λ, ∂hη=∂hη)=λ.*x.-A'∂hη.(y.-A*x)

tₖ(k:: Int64)=1/β #Stepsize rule p/ subgradient descent

include("../Métodos/subgradient_descent_plot.jl")

p₁=subgradient_descent(f, ∇f, tₖ, copy(x₀), k_max, ϵ)

include("../Métodos/acc_subgradient_descent_plot.jl")

p₂=acc_subgradient_descent(f, ∇f, β, copy(x₀), k_max, ϵ)

plot(p₁, p₂)