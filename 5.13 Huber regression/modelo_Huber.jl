include("variables_Huber.jl")

hη(w:: Float64; η=η)=(abs(w)>η)*η*(abs(w)-η/2)+!(abs(w)>η)*w^2/2
f(x:: Vector{Float64}; A=A, y=y, λ=λ, hη=hη)=sum(hη.(y.-A*x))+λ*norm(x)^2/2

∂hη(w:: Float64; η=η)=(abs(w)>η)*η*sign(w)+!(abs(w)>η)*w
∇f(x:: Vector{Float64}; A=A, y=y, λ=λ, ∂hη=∂hη)=λ.*x.-A'∂hη.(y.-A*x)

include("../Métodos/gradient_descent_plot.jl")

p₁=gradient_descent(f, ∇f, β, x₀, t_max)

include("../Métodos/acc_gradient_descent_plot.jl")

p₂=acc_gradient_descent(f, ∇f, β, x₀, t_max)

plot(p₁, p₂)