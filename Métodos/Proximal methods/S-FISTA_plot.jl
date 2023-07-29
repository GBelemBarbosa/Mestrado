using Plots
using LaTeXStrings

function SFISTA(f:: Function, g:: Function, h:: Function, hμ:: Function, ∇Fμ:: Function, prox_gL:: Function, Lf:: Number, μ:: Number, α:: Number, x₀:: Array{<:Number}, k_max:: Int64, ϵ:: Number) 
    y, x=x₀, x₀
    t=1
    fx=f(x)+g(x)
    L̃=Lf+α/μ
    hist=[fx+h(x), fx+hμ(x)]
    
    for k=0:k_max
        ∇Fμy=∇Fμ(y)

        if norm(∇Fμy, Inf)<ϵ
            break
        end

        x_, x=x, prox_gL(L̃, y.-∇Fμy./L̃)
        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_)

        fx=f(x)+g(x)
        push!(hist, fx+h(x))
        push!(hist, fx+hμ(x))
    end 

    fx=f(x)+g(x)
    println(norm(∇Fμ(y), Inf), " ", fx+h(x), " ", fx+hμ(x))
    x, scatter(eachindex(hist[begin:2:end]), [hist[begin:2:end], hist[2:2:end]], 
                title=L"F(x^{(k)})"*" e "*L"F_\mu(x^{(k)})",
                label=[L"F(x^{(k)})" L"F_\mu(x^{(k)})"])
end