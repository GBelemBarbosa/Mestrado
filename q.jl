using Plots, LaTeXStrings

A=[4 -2 1
1 -1 1
0 0 1
1 1 1
4 2 1]

y₀=[10/18
10/13
1
10/15
10/16]

y=log.(y₀)

α=(A'A)\(A'y)

b=1/sqrt(-α[1])
a=-α[2]*b^2/2

x=[-2
-1 
0 
1 
2]

p=plot(x, y₀, label=L"y_i")

ϕ₂(x:: Number)=exp(-((a-x)/b)^2)

plot!(x, ϕ₂, label=L"aproximação")

aux(x:: Number)=log(ϕ₂(x))

q=plot(x, y, label=L"y_i")

plot!(x, aux, label=L"aproximação")