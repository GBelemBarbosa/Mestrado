using Plots
using StatsPlots
using LaTeXStrings
using LinearAlgebra

include("parameters.jl")
include("functions.jl")
include("first_order_generalized_descent.jl")
include("second_order_generalized_descent.jl")
include("newton.jl")

nam=repeat(vcat("0.25", ["2^-"*string(i) for i=12:10:32]), outer = 4)
ctg=repeat(["Gradiente", "Newton", "CP1", "DFP"], inner=length(ϵ_vec))

convergence=[reduce(vcat, [[(j==1)*n_vec[k]+(j==2)*norm_vec[i]+0.0 for i=eachindex(norm_vec), j=1:6] for k=eachindex(n_vec)]) for p=eachindex(f_vec)]

for j=eachindex(f_vec)
    global f=f_vec[j]
    global ∇f=∇_vec[j]
    global ∇²f=∇²_vec[j]
    xₒₚₜ_at=xₒₚₜ_vec[j]
    fxₒₚₜ_at=fxₒₚₜ_vec[j]
    println(string(f))

    for k=eachindex(n_vec)
        global n=n_vec[k]
        if j==2 && n>10
            break
        end
        println("n=",n)

        for p=eachindex(norm_vec)
            r=norm_vec[p]
            println("r=", r)
            distxₒₚₜ_grad, distxₒₚₜ_newton, distxₒₚₜ_CP1, distxₒₚₜ_DFP=[zeros(k_max) for i=1:4]
            distx₀=0
            distfxₒₚₜ_grad, distfxₒₚₜ_newton, distfxₒₚₜ_CP1, distfxₒₚₜ_DFP=[zeros(k_max) for i=1:4]
            distfx₀=0
            tk_grad, tk_newton, tk_CP1, tk_DFP=[zeros(k_max) for i=1:4]
            converge_grad, converge_newton, converge_CP1, converge_DFP=[[0.0] for i=1:4]
            kϵ_grad, kϵ_newton, kϵ_CP1, kϵ_DFP=[zeros(length(ϵ_vec)) for i=1:4]
            global xₒₚₜ=xₒₚₜ_at()
            global fxₒₚₜ=fxₒₚₜ_at()
            global H₀=diagm([1 for i=1:n])

            for i=1:s_max
                x₀=randn(n).*r
                distx₀+=norm(x₀.-xₒₚₜ, p)/s_max
                distfx₀+=(f(x₀)-fxₒₚₜ)/s_max
                
                print(i)
                FOGD(x₀, distxₒₚₜ_grad, distfxₒₚₜ_grad, tk_grad, converge_grad, kϵ_grad)
                newton(x₀, distxₒₚₜ_newton, distfxₒₚₜ_newton, tk_newton, converge_newton, kϵ_newton)
                SOGD(Hₖ_CP1, x₀, distxₒₚₜ_CP1, distfxₒₚₜ_CP1, tk_CP1, converge_CP1, kϵ_CP1)
                SOGD(Hₖ_DFP, x₀, distxₒₚₜ_DFP, distfxₒₚₜ_DFP, tk_DFP, converge_DFP, kϵ_DFP)
                print("\e[2K") # clear whole line
                print("\e[1G") # move cursor to column 1   
            end

            println("% Newton: ", converge_newton[1])
            println("% CP1: ", converge_CP1[1])
            println("% DFP: ", converge_DFP[1])
            convergence[j][(k-1)*3+p, 3:6]=[converge_grad[1], converge_newton[1], converge_CP1[1], converge_DFP[1]]

            pushfirst!(distxₒₚₜ_grad, distx₀)
            pushfirst!(distxₒₚₜ_newton, distx₀)
            pushfirst!(distxₒₚₜ_CP1, distx₀)
            pushfirst!(distxₒₚₜ_DFP, distx₀)

            pushfirst!(distfxₒₚₜ_grad, distfx₀)
            pushfirst!(distfxₒₚₜ_newton, distfx₀)
            pushfirst!(distfxₒₚₜ_CP1, distfx₀)
            pushfirst!(distfxₒₚₜ_DFP, distfx₀)

            plot_max=Int64(round(max(kϵ_grad[1], kϵ_newton[1], kϵ_CP1[1], kϵ_DFP[1])))

            c=min(1, 1/log(distx₀))
            distxₒₚₜ_plot=plot(0:plot_max-1, distxₒₚₜ_CP1[1:plot_max].^c, 
                label="CP1",
                title=L"||\overline{x^k}-x^*||_\infty^c,\ r="*string(r)*", n="*string(n),
                xlabel=L"k",
                ylabel=L"||\overline{x^k}-x^*||_\infty^c")
            plot!(0:plot_max-1, distxₒₚₜ_DFP[1:plot_max].^c, 
                label="DFP")
            plot!(0:plot_max-1, distxₒₚₜ_grad[1:plot_max].^c, 
                label="Gradiente")
            plot!(0:plot_max-1, distxₒₚₜ_newton[1:plot_max].^c, 
                label="Newton")
            
            c=min(1, 1/log(distfx₀))
            distfxₒₚₜ_plot=plot(0:plot_max-1, abs.(distfxₒₚₜ_CP1[1:plot_max]).^c, 
                label="CP1",
                title=L"(f(\overline{x^k})-f(x^*))^c,\ r="*string(r)*", n="*string(n),
                xlabel=L"k",
                ylabel=L"(f(\overline{x^k})-f(x^*))^c")
            plot!(0:plot_max-1, abs.(distfxₒₚₜ_DFP[1:plot_max]).^c, 
                label="DFP")
            plot!(0:plot_max-1, abs.(distfxₒₚₜ_grad[1:plot_max]).^c, 
                label="Gradiente")
            plot!(0:plot_max-1, abs.(distfxₒₚₜ_newton[1:plot_max]).^c, 
                label="Newton")

            tk_plot=plot(1:plot_max, tk_CP1[1:plot_max], 
                label="CP1",
                title=L"\overline{t_k},\ r="*string(r)*", n="*string(n),
                xlabel=L"k",
                ylabel=L"\overline{t_k}")
            plot!(1:plot_max, tk_DFP[1:plot_max], 
                label="DFP")
            plot!(1:plot_max, tk_grad[1:plot_max], 
                label="Gradiente")
            plot!(1:plot_max, tk_newton[1:plot_max], 
                label="Newton")
            
            kϵ_plot=groupedbar(nam, hcat(kϵ_grad, kϵ_newton, kϵ_CP1, kϵ_DFP),
                title=L"k"*" para obter "*L"||\nabla f(x^k)||_\infty <\epsilon_i,\ r="*string(r)*", n="*string(n),
                group=ctg,
                xlabel=L"\epsilon_i",
                ylabel=L"k")

            png(distxₒₚₜ_plot, "distx_"*string(n)*"_"*string(f)*string(r))
            png(distfxₒₚₜ_plot, "distf_"*string(n)*"_"*string(f)*string(r))
            png(tk_plot, "tk_"*string(n)*"_"*string(f)*string(r))
            png(kϵ_plot, "k_"*string(n)*"_"*string(f)*string(r))
        end
    end
end

for i=1:4
    println("\\begin{table}[H]
    \\centering
    \\begin{tabular}{@{}lccccc}
    \\hline
    \$n\$& \$r\$& Gradiente& Newton& CP1& DFP\\\\
    \\hline")
    for j=1:size(convergence[i])[1]
        if i!=2 || convergence[i][j, 1]<=10
            for k=1:2
                print(string(Int64(convergence[i][j, k]))*"& ")
            end
            for k=3:5
                print(string(round(convergence[i][j, k]; digits=2))*"& ")
            end
            println(string(round(convergence[i][j, 6]; digits=2))*"\\\\")
        end
    end
    println("\\hline
    \\end{tabular}
    \\caption{Proporção de obtenção da solução ótima para \$f\$ "*string(f_vec[i])*".}
    \\end{table}")
    println()
end