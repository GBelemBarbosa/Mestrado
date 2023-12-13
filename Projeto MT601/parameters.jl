n_vec=[2, 10, 100]
ϵ_vec=[2.0^-i for i=2:10:32]
norm_vec=[1, 10, 100]
s_max=100
p=Inf

k_max=10000
α=10^-4
β=10^-3
γ=10^-6
σ=1/2
ρ=10^-3

t₀=1
tₖ(k:: Int64, t:: Number, x:: Array{<:Number}, d:: Array{<:Number}; σ=σ)=σ*t
resize(k:: Int64, t_:: Number, t:: Number)=1
dₖ(k:: Int64, x_:: Array{<:Number}, x:: Array{<:Number}, ∇fx_:: Array{<:Number}, ∇fx:: Array{<:Number}, d:: Array{<:Number})=.-∇fx

function Hₖ_CP1(x_:: Array{<:Number}, x:: Array{<:Number}, ∇fx_:: Array{<:Number}, ∇fx:: Array{<:Number}, H_:: Array{<:Number}, H:: Array{<:Number})
    yₖ=∇fx.-∇fx_
    wₖ=x.-x_.-H*yₖ
    wₖTyₖ=wₖ'yₖ

    if wₖTyₖ>=0
        return H.+wₖ*(wₖ'./wₖTyₖ)
    end

    return H
end

function Hₖ_DFP(x_:: Array{<:Number}, x:: Array{<:Number}, ∇fx_:: Array{<:Number}, ∇fx:: Array{<:Number}, H_:: Array{<:Number}, H:: Array{<:Number})
    sₖ=x.-x_
    yₖ=∇fx.-∇fx_
    sₖTyₖ=sₖ'yₖ

    if sₖTyₖ>=0
        zₖ=H*yₖ
        return H.+sₖ*(sₖ'./sₖTyₖ).-zₖ*(zₖ'./zₖ'yₖ)
    end

    return H
end

fail(H_:: Array{<:Number}, H:: Array{<:Number})=diagm(ones(size(H)[1]))