include("../Métodos/Proximal methods hist/BB.jl")
include("../Métodos/Proximal methods hist/nmBB.jl")
include("../Métodos/Proximal methods hist/FISTA.jl")
include("../Métodos/Proximal methods hist/proximal_subgradient.jl")
include("../Métodos/Proximal methods hist/mAPG.jl")
include("../Métodos/Proximal methods hist/nmAPG.jl")
include("../Métodos/Proximal methods hist/mAPGLS.jl")
include("../Métodos/Proximal methods hist/nmAPGLS.jl")
include("../Métodos/Proximal methods hist/nmAPGLS2.jl")
include("../Métodos/Proximal methods hist/IBPGLS.jl")
include("../Métodos/Proximal methods hist/newAPG_vs.jl")
include("../Métodos/Proximal methods hist/newAPG_vs_x_2.jl")
include("../Métodos/Proximal methods hist/newAPG_vs_y_2.jl")

const α₀ = norm(∇f(x₀).-∇f(x₀.+10^-5))/(sqrt(n)*10^-5)
x_BB, histF, histnψ = BB(F, ∇f, pαₖ, x₀, α₀, 10^-30, 10^30, 5, k_max; ϵ=ϵ)   
plot!(pltF, histF, marker=:circle, label="BB")
plot!(pltnψ, histnψ, marker=:circle, label="BB")
println("x_BB: ", 1-norm(x_BB, 0)/n, " & ", f(x_BB))

x_mBB, histF, histnψ, ls_mBB = nmBB(F, ∇f, pαₖ, x₀, 4, 0.01, 10^-30, 10^30, 1, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:circle, label="mBB")
plot!(pltnψ, histnψ, marker=:circle, label="mBB")
println("x_mBB: ", 1-norm(x_mBB, 0)/n, " & ", f(x_mBB))

x_nmBB, histF, histnψ, ls_nmBB = nmBB(F, ∇f, pαₖ, x₀, 4, 0.01, 10^-30, 10^30, 5, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:circle, label="nmBB")
plot!(pltnψ, histnψ, marker=:circle, label="nmBB")
println("x_nmBB: ", 1-norm(x_nmBB, 0)/n, " & ", f(x_nmBB))

x_mAPG, histF, histnψ = mAPG(F, ∇f, proxα, x₀, 1/L, 1/L, k_max; ϵ=ϵ) # 2 prox steps per iteration
plot!(pltF, histF, marker=:diamond, label="mAPG")
plot!(pltnψ, histnψ, marker=:diamond, label="mAPG")
println("x_mAPG: ", 1-norm(x_mAPG, 0)/n, " & ", f(x_mAPG))

x_nmAPG, histF, histnψ, ls_nmAPG = nmAPG(F, ∇f, proxα, x₀, 1/L, 1/L, 0.8, 10^-4, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:diamond, label="nmAPG")
plot!(pltnψ, histnψ, marker=:diamond, label="nmAPG")
println("x_nmAPG: ", 1-norm(x_nmAPG, 0)/n, " & ", f(x_nmAPG))

x_mAPGLS, histF, histnψ, ls_mAPGLS = mAPGLS(F, ∇f, proxα, x₀, 2/5, 10^-4, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:diamond, label="mAPGLS")
plot!(pltnψ, histnψ, marker=:diamond, label="mAPGLS")
println("x_mAPGLS: ", 1-norm(x_mAPGLS, 0)/n, " & ", f(x_mAPGLS))

x_nmAPGLS, histF, histnψ, ls_nmAPGLS = nmAPGLS(F, ∇f, proxα, x₀, 2/5, 0.8, 10^-4, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:diamond, label="nmAPGLS")
plot!(pltnψ, histnψ, marker=:diamond, label="nmAPGLS")
println("x_nmAPGLS: ", 1-norm(x_nmAPGLS, 0)/n, " & ", f(x_nmAPGLS))

x_nmAPGLS2, histF, histnψ, ls_nmAPGLS2 = nmAPGLS2(F, ∇f, proxα, x₀, 2/5, 0.8, 10^-4, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:diamond, label="nmAPGLS2")
plot!(pltnψ, histnψ, marker=:diamond, label="nmAPGLS2")
println("x_nmAPGLS2: ", 1-norm(x_nmAPGLS2, 0)/n, " & ", f(x_nmAPGLS2))

x_PG, histF, histnψ = PG(F, ∇f, Lₖ, x₀, L, k_max; ϵ=ϵ)
plot!(pltF, histF, linestyle=:dash, label="PG")
plot!(pltnψ, histnψ, linestyle=:dash, label="PG")
println("x_PG: ", 1-norm(x_PG, 0)/n, " & ", f(x_PG))

x_FISTA, histF, histnψ = FISTA(F, ∇f, Lₖ, x₀, L, k_max; ϵ=ϵ)
plot!(pltF, histF, linestyle=:dash, label="FISTA")
plot!(pltnψ, histnψ, linestyle=:dash, label="FISTA")
println("x_FISTA: ", 1-norm(x_FISTA, 0)/n, " & ", f(x_FISTA))

α°ₖ(k:: Int64, αₖ:: Number)=9
β°ₖ(k:: Int64, βₖ:: Number)=0.9
λ°ₖ(k:: Int64, λₖ:: Number)=0.5 
const λₘᵢₙ=0.98*(1-0.1)/Lf
const d=((1-0.1)/λₘᵢₙ-Lf)/2
pₖ(k:: Int64)=1
#=x_IBPGmLS, histF, histnψ, ls_IBPGmLS = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:utriangle, label="IBPGmLS")
plot!(pltnψ, histnψ, marker=:utriangle, label="IBPGmLS")
println("x_IBPGmLS: ", 1-norm(x_IBPGmLS, 0)/n, " & ", f(x_IBPGmLS))=#

λ°ₖLf(k:: Int64, λₖ:: Number; Lf=Lf)=0.99*(1-0.1)/Lf
x_IBPGmLSLf, histF, histnψ, ls_IBPGmLSLf = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖLf, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:utriangle, label="IBPGmLSLf")
plot!(pltnψ, histnψ, marker=:utriangle, label="IBPGmLSLf")
println("x_IBPGmLSLf: ", 1-norm(x_IBPGmLSLf, 0)/n, " & ", f(x_IBPGmLSLf))

pₖ(k:: Int64)=0.7
#=x_IBPGnmLS, histF, histnψ, ls_IBPGnmLS = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖ, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:utriangle, label="IBPGnmLS")
plot!(pltnψ, histnψ, marker=:utriangle, label="IBPGnmLS")
println("x_IBPGnmLS: ", 1-norm(x_IBPGnmLS, 0)/n, " & ", f(x_IBPGnmLS))=#

x_IBPGnmLSLf, histF, histnψ, ls_IBPGnmLSLf = IBPGLS(F, ∇f, ℘hλg, α°ₖ, β°ₖ, λ°ₖLf, pₖ, x₀, 0.4, 0.35, 0.45, d, 0.1, λₘᵢₙ, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:utriangle, label="IBPGnmLSLf")
plot!(pltnψ, histnψ, marker=:utriangle, label="IBPGnmLSLf")
println("x_IBPGnmLSLf: ", 1-norm(x_IBPGnmLSLf, 0)/n, " & ", f(x_IBPGnmLSLf))

function γ(k:: Int64, t_:: Number)
    t=(1+sqrt(1+4*t_^2))/2

    return t, (t_-1)/t
end
Q(k:: Int64)=0.99^k
E(k:: Int64)=k^-1.1
const λ₁=1/α₀
x_newAPG_vs, histF, histnψ, ls_newAPG_vs = newAPG_vs(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:x, label="newAPG_vs")
plot!(pltnψ, histnψ, marker=:x, label="newAPG_vs")
println("x_newAPG_vs: ", 1-norm(x_newAPG_vs, 0)/n, " & ", f(x_newAPG_vs))

x_newAPG_vs_x_2, histF, histnψ, ls_newAPG_vs_x_2 = newAPG_vs_x_2(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:x, label="newAPG_vs_x_2")
plot!(pltnψ, histnψ, marker=:x, label="newAPG_vs_x_2")
println("x_newAPG_vs_x_2: ", 1-norm(x_newAPG_vs_x_2, 0)/n, " & ", f(x_newAPG_vs_x_2))

x_newAPG_vs_y_2, histF, histnψ, ls_newAPG_vs_y_2 = newAPG_vs_y_2(F, h, ∇f, Tλ, γ, Q, E, x₀, λ₁, 0.99, 0.95, 10^4, 0.8, k_max; ϵ=ϵ)
plot!(pltF, histF, marker=:x, label="newAPG_vs_y_2")
plot!(pltnψ, histnψ, marker=:x, label="newAPG_vs_y_2")
println("x_newAPG_vs_y_2: ", 1-norm(x_newAPG_vs_y_2, 0)/n, " & ", f(x_newAPG_vs_y_2))

pltF