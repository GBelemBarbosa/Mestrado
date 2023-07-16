using Plots
using LaTeXStrings

function SFISTA(f:: Function, g:: Function, h:: Function, hμ:: Function, ∇f:: Function, prox_gL:: Function, Lf:: Number, μ:: Number, α:: Number, x₀:: Array{T, N}, k_max:: Int64, ϵ:: Number) where {T, N}
    y, x=x₀, x₀
    t=1
    fx=f(x)+g(x)
    L=Lf+α/μ
    hist=[fx+h(x), fx+hμ(x)]
    
    for k=0:k_max
        ∇fx=∇f(x)

        if norm(∇fx, Inf)<ϵ
            break
        end

        x_, x=x, prox_gL(x.-∇fx./L)
        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_)

        fx=f(x)+g(x)
        push!(hist, fx+h(x))
        push!(hist, fx+hμ(x))
    end 

    fx=f(x)+g(x)
    println(norm(∇f(x), Inf), " ", fx+h(x), " ", fx+hμ(x))
    scatter(eachindex(hist[begin:2:end]), [hist[begin:2:end], hist[2:2:end]], 
                title=L"F(x^{(k)})"*" e "*L"F_\mu(x^{(k)})",
                label=[L"F(x^{(k)})" L"F_\mu(x^{(k)})"])
end