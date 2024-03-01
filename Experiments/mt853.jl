"""
Simple and naive implementations of optimization methods for the Optimization
for Data Science course.
"""

using MAT
using LinearAlgebra
using SparseArrays

"""
    struct LinearL1Problem

    struct that groups the information of a linear problem with ℓ_1 regularization.

        - A: matrix in CSC format;
        - b: right-hand side of least squares problem;
        - λ: wright to balance fitting with sparsity of the solution;
        - ftarget: A target for the Lasso objetive tha ensure a good enough solution.
"""
struct LinearL1Problem
    A::SparseMatrixCSC{Float64}
    b::Vector{Float64}
    λ::Float64
    ftarget::Float64
end

"""
    readlasso(filename)

    Read a lasso problem from file.
"""
function readlasso(filename)
    vars = matread(filename)
    return LinearL1Problem(vars["A"], vec(vars["b"]), vars["lambda"], vars["ftarget"])
end

"""
    readlogreg(filename)

    Read a ℓ1 regularized logistic regression problem from file.
"""
function readlogreg(filename)
    vars = matread(filename)
    return LinearL1Problem(vars["A"], vec(vars["b"]), vars["lambdalog"], vars["flogtarget"])
end

"""
    grad_descent(f, gradf, stepsize; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)

    Gradient descent method. 

        - x0: starting point;
        - f: objective function;
        - gradf: function to compute objective and its gradient;
        - stepsize: function to compute the stepsize;
        - ϵ: method stop if | ∇f(x) | <= ϵ;
        - f_target: stop if f(x) <= f_target;
        - max_iter: maximum number of iterations
"""
function grad_descent(x0, f, gradf, stepsize; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)
    x = copy(x0)
    fval, ∇f = gradf(x)
    histf = [fval]
    hist∇f = [norm(∇f, Inf)]
    iter = 0
    while hist∇f[end] > ϵ && fval > ftarget && iter < max_iter
        d = -∇f
        α = stepsize(x, d, fval, ∇f, f, gradf)
        @. x = x + α * d
        fval, ∇f = gradf(x)
        iter += 1
        append!(histf, fval)
        append!(hist∇f, norm(∇f, Inf))
    end
    return x, fval, ∇f, histf, hist∇f
end

"""
    spectral_grad(f, gradf; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)

    Spectral Gradient method for quadratics (there is no globalization)

        - x0: starting point;
        - f: objective function;
        - gradf: function to compute objective and its gradient;
        - ϵ: method stop if | ∇f(x) | <= ϵ;
        - f_target: stop if f(x) <= f_target.
        - max_iter: maximum number of iterations
"""
function spectral_grad(x0, f, gradf; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)
    x = copy(x0)
    fval, ∇f = gradf(x)
    histf = [fval]
    hist∇f = [norm(∇f, Inf)]

    lastx = copy(x)
    lastg = copy(∇f)
    Δx = x - lastx
    Δg = ∇f - lastg
    iter = 0
    while hist∇f[end] > ϵ && fval > ftarget && iter < max_iter
        if iter == 0
            α = 1.0e-5 * norm(∇f)
        else
            @. Δx = x - lastx
            @. Δg = ∇f - lastg
            α = dot(Δx, Δx) / dot(Δx, Δg)
        end
        lastx .= x
        lastg .= ∇f
        @. x = x - α * ∇f
        fval, ∇f = gradf(x)
        iter += 1
        append!(histf, fval)
        append!(hist∇f, norm(∇f, Inf))
    end
    return x, fval, ∇f, histf, hist∇f
end


"""
    proxgrad(f, gradf, stepsize; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)

    Proximal gradient descent method to minimize

    f(x) + h(x)

        - x0: starting point;
        - f: first term of the objective function;
        - gradf: function to compute objective and its gradient;
        - h: second term of the objetive;
        - prox!: proximal operator associated with function h;
        - α: the cte stepsize;
        - f_target: stop if f(x) <= f_target;
        - max_iter: maximum number of iterations
"""
function proxgrad(x0, f, gradf, h, prox!, α; ftarget=-1.0e20, max_iter=2000, adapt=false)
    x = copy(x0)
    fval, ∇f = gradf(x)
    fval += h(x)
    histf = [fval]
    hist∇f = [norm(∇f, Inf)]
    iter = 0
    while fval > ftarget && iter < max_iter
        d = -∇f
        lastx = copy(x)
        lastfval = fval
        @. x = x + α * d
        prox!(x, α)
        fval, ∇f = gradf(x)
        fval += h(x)
        if adapt && fval > lastfval - 1.0 / (2.0 * α) * norm(x - lastx)^2
            # println("Decreasing α.")
            α /= 10
            x = lastx
            fval, ∇f = gradf(x)
            fval += h(x)
        end
        iter += 1

        append!(histf, fval)
        append!(hist∇f, norm(∇f, Inf))
    end
    return x, fval, ∇f, histf, hist∇f
end


"""
    spectral_grad(f, gradf; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)

    Spectral proximal gradient method (there is no globalization) to minimize

    f(x) + h(x)

        - x0: starting point;
        - f: objective function;
        - gradf: function to compute objective and its gradient;
        - h: second term of the objetive;
        - prox!: proximal operator associated with function h;
        - f_target: stop if f(x) <= f_target.
        - max_iter: maximum number of iterations
"""
function spectral_proxgrad(x0, f, gradf, h, prox!; ftarget=-1.0e20, max_iter=2000)
    x = copy(x0)
    fval, ∇f = gradf(x)
    fval += h(x)
    histf = [fval]
    hist∇f = [norm(∇f, Inf)]

    m = 10
    fmem = Inf * ones(m)
    fmem[1] = fval
    lastx = copy(x) .+ 1
    lastg = copy(∇f)
    Δx = x - lastx
    Δg = ∇f - lastg
    iter = 0
    while norm(Δx) > 0.0 && fval > ftarget && iter < max_iter
        if iter == 0
            α = 1.0e-5 / norm(∇f)
        else
            α = dot(Δx, Δx) / dot(Δx, Δg)
        end
        lastx .= x
        lastg .= ∇f
        fmax = maximum(fmem)
        fval = Inf
        ntries = 0
        while (fval >= fmax)
            @. x = lastx - α * lastg
            prox!(x, α)
            fval, ∇f = gradf(x)
            fval += h(x)
            α /= 10.0
            ntries += 1
        end
        iter += 1
        fmem[iter%m+1] = fval
        for i = 1:ntries
            append!(histf, fval)
            append!(hist∇f, norm(∇f, Inf))
        end
        @. Δx = x - lastx
        @. Δg = ∇f - lastg
    end
    return x, fval, ∇f, histf, hist∇f
end


"""
    heavy_ball(f, gradf, L, σ; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)

    Heavy ball method. 

        - x0: starting point;
        - f: objective function;
        - gradf: function to compute objective and its gradient;
        - L: L-smooth cte;
        - σ: strong convexity cte;
        - ϵ: method stop if | ∇f(x) | <= ϵ;
        - f_target: stop if f(x) <= f_target;
        - max_iter: maximum number of iterations
"""
function heavy_ball(x0, f, gradf, L, σ; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)
    α = 4.0 / (sqrt(L) + sqrt(σ))^2
    θ = max((1 - sqrt(α * L))^2, (1 - sqrt(α * σ))^2)
    # Only the first step is 1/L to ensure decrease before having momentum.
    α = 1 / L
    x = copy(x0)
    fval, ∇f = gradf(x)
    histf = [fval]
    hist∇f = [norm(∇f, Inf)]
    v = zeros(length(x))
    iter = 0
    while hist∇f[end] > ϵ && fval > ftarget && iter < max_iter
        v = θ * v - α * ∇f
        @. x = x + v
        α = 4.0 / (sqrt(L) + sqrt(σ))^2
        fval, ∇f = gradf(x)
        iter += 1
        append!(histf, fval)
        append!(hist∇f, norm(∇f, Inf))
    end
    return x, fval, ∇f, histf, hist∇f
end


"""
    optimal_method(f, gradf, L; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000)

    Optimal method. 

        - x0: starting point;
        - f: objective function;
        - gradf: function to compute objective and its gradient;
        - L: L-smooth cte;
        - ϵ: method stop if | ∇f(x) | <= ϵ;
        - f_target: stop if f(x) <= f_target;
        - max_iter: maximum number of iterations
"""
function optimal_method(x0, f, gradf, L; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000, restart=false)
    α = 1 / L
    x = copy(x0)
    xprev = copy(x)
    fval, ∇f = gradf(x)
    histf = [fval]
    hist∇f = [norm(∇f, Inf)]
    v = zeros(length(x))
    iter = 0
    momentum = iter
    while hist∇f[end] > ϵ && fval > ftarget && iter < max_iter
        if momentum == 0
            θ = 0
        else
            θ = (momentum - 1) / (momentum + 2)
        end
        y = x + θ * (x - xprev)
        fvaly, ∇fy = gradf(y)
        @. xprev = x
        @. x = y - α * ∇fy
        fval, ∇f = gradf(x)
        iter += 1
        append!(histf, fval)
        append!(hist∇f, norm(∇f, Inf))
        if restart && histf[end] > histf[end-1]
            @. xprev = x
            momentum = 0
        else
            momentum += 1
        end
    end
    return x, fval, ∇f, histf, hist∇f
end


"""
    naive_incremental_grad(x0, gradf, m; epochs=100)

    I am going to be lazy and use stepsize = 1/L. 
"""
function naive_incremental_grad(x0, gradf, m, L; epochs=100)
    x = copy(x0)
    for epoch = 1:epochs
        window = 10
        α = 1 / L #(L*div(epoch + window - 1, window))
        for i = 1:m
            g = gradf(x, i)
            @. x = x - α * g
        end
    end
    return x
end


function naive_fista(x0, f, gradf, h, prox!, L; ϵ=1.0e-5, ftarget=-1.0e20, max_iter=2000, restart=false)
    α = 1 / L
    x = copy(x0)
    y = copy(x)
    xprev = copy(x)
    fval, ∇f = gradf(y)
    histf = [fval + h(y)]
    hist∇f = [norm(∇f, Inf)]
    iter = 0
    θ = 1
    while hist∇f[end] > ϵ && histf[end] > ftarget && iter < max_iter
        @. xprev = x
        @. x = y - α * ∇f
        prox!(x, α)
        θprev = θ
        θ = 0.5 * (1 + sqrt(1 + 4 * θ^2))
        @. y = x + (θprev - 1) / θ * (x - xprev)
        fval, ∇f = gradf(y)
        iter += 1
        append!(histf, fval + h(y))
        append!(hist∇f, norm(∇f, Inf))
        if restart && histf[end] > histf[end-1]
            @. xprev = x
            @. y = x
            θ = 1
        end
    end
    return y, histf[end], ∇f, histf, hist∇f
end

function scalar_soft_thresholding(x, α)
    if x < -α
        return x + α
    elseif x <= α
        return 0.0
    else
        return x - α
    end
end


function mydot(A, i, res)
    dotval = 0.0
    @inbounds for j = A.colptr[i]:A.colptr[i+1]-1
        dotval += A.nzval[j] * res[A.rowval[j]]
    end
    return dotval
end


function update_res!(res, A, i, delta)
    @inbounds for j = A.colptr[i]:A.colptr[i+1]-1
        res[A.rowval[j]] += delta * A.nzval[j]
    end
end



function lasso_coord(A, b, λ, ftarget; epochs=100)
    m, n = size(A)

    # Compute L_i
    L = Vector{Float64}(undef, n)
    for i = 1:n
        L[i] = norm(A[:, i])^2
    end

    # Start from 0 so the initial residual is easy
    x = zeros(n)
    res = -b

    epoch = 0
    f = 0.5 * norm(b)^2
    histf = [f]
    @inbounds while epoch < epochs && f > ftarget
        for i = 1:n
            new_coord = x[i] - 1 / L[i] * mydot(A, i, res)
            new_coord = scalar_soft_thresholding(new_coord, λ / L[i])
            delta = new_coord - x[i]
            update_res!(res, A, i, delta)
            x[i] = new_coord
        end

        # Update function value
        f = 0.5 * norm(res)^2 + λ * norm(x, 1)
        epoch += 1
        append!(histf, f)
    end
    return x, f, histf
end







