using Plots
using LaTeXStrings

function SFISTA(f:: Function, g:: Function, h:: Function, hμ:: Function, ∇Fμ:: Function, prox_gL:: Function, Lf:: Number, μ:: Number, α:: Number, x₀:: Array{<:Number}, k_max:: Int64; ϵ=eps(), p=Inf) 
    y=x_=x=x₀
    t=1
    fx=f(x)+g(x)
    L̃=Lf+α/μ
    hist=[fx+h(x), fx+hμ(x)]
    
    k=1
    while true
        ∇Fμy=∇Fμ(y)
        x=prox_gL(L̃, y.-∇Fμy./L̃)

        fx=f(x)+g(x)
        push!(hist, fx+h(x))
        push!(hist, fx+hμ(x))

        if norm(∇Fμ(x).-∇Fμy.+(y.-x).*L̃, p)<ϵ || k==k_max
            break
        end
        k+=1

        t_, t=t, (1+sqrt(1+4*t^2))/2
        y=x.+((t_-1)/t).*(x.-x_)
        x_=x
    end 

    println(norm(∇Fμ(y), p), " ", hist[end-1], " ", hist[end])
    x, scatter(eachindex(hist[begin:2:end]), [hist[begin:2:end], hist[2:2:end]], 
                title=L"F(x^{(k)})"*" e "*L"F_\mu(x^{(k)})",
                label=[L"F(x^{(k)})" L"F_\mu(x^{(k)})"])
end