using LinearAlgebra, Plots, LaTeXStrings

A=[-2.12 1.1 0 0 0 
 0.9 -2.12 1.1 0 0 
 0 0.9 -2.12 1.1 0 
 0 0 0.9 -2.12 1.1 
 0 0 0 2 -3]

b=vcat([-0.04*i for i=0.2:0.2:0.8], 0.4)

y=A\b

pushfirst!(y, 0)

f(x:: Number)=((3 + sqrt(13))*exp(sqrt(13))*(3*x + 1) + (sqrt(13) - 3)*(3*x + 1) + 40*exp(- (1 + sqrt(13))*(x - 1)/2) - (sqrt(13) - 3)*exp((sqrt(13) - 1)*x/2) - 40*exp(((sqrt(13) - 1)*x + sqrt(13) + 1)/2) - (3 + sqrt(13))*exp(sqrt(13) - (1 + sqrt(13))*x/2))/(9*(-3 + sqrt(13) + (3 + sqrt(13))*exp(sqrt(13))))

p=plot(0:0.2:1, y, label=L"y_i"; marker=(:circle,5))

plot!(0:0.01:1, f, label=L"y(x)")

png(p, "plot")