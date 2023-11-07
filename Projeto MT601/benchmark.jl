using Plots
using Manopt, Manifolds
using LinearAlgebra
using BenchmarkTools

include("parameters.jl")
include("functions_manifold.jl")

gradient(x:: Array{<:Number}; f=f, ∇f=∇f, n=n)=gradient_descent(Euclidean(n), f, ∇f, x; stopping_criterion=StopWhenAny(StopAfterIteration(1000), StopWhenGradientNormLess(10^(-6))))
CP1(x:: Array{<:Number}; f=f, ∇f=∇f, n=n)=quasi_Newton(Euclidean(n), f, ∇f, x; direction_update=InverseSR1())
DFP(x:: Array{<:Number}; f=f, ∇f=∇f, n=n)=quasi_Newton(Euclidean(n), f, ∇f, x; direction_update=InverseDFP())

function repeval(method:: Function; f=f, n=n, s_max=s_max)
    for i=1:s_max
        method(randn(n).*10)
    end
end

for j=eachindex(f_vec)
    global f=f_vec[j]
    global ∇f=∇_vec[j]
    global ∇²f=∇²_vec[j]

    println(string(f))

    for k=eachindex(n_vec)
        global n=n_vec[k]
        
        println("n=", n)

        println("Gradient:")
        @btime repeval(gradient)

        println("CP1:")
        @btime repeval(CP1)

        println("DFP:")
        @btime repeval(DFP)
    end

    println()    
end