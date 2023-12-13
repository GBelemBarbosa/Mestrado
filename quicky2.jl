using LinearAlgebra, Plots

b=[24.9; 25; 25.1; 21.5; 19.5; 18.3; 19.4; 22.8; 24.5; 25.3]

ω=pi/6
A=reduce(vcat, [[1, cos(ω*t), sin(ω*t)]' for t=1:10])

α=(A'A)\(A'b) #Solução do sistema normal

T(t:: Number; α=α, ω=ω)=α[1]+α[2]*cos(ω*t)+α[3]*sin(ω*t) #Definição da função da curva

res=sum((T.(1:10).-b).^2) #Resíduo

T_nov=T(11) #Temperatura média em novembro
T_dez=T(12) #Temperatura média em dezembro

p=plot(1:0.1:12, T, label="Curva ajustada", legend=:bottomright)
scatter!(1:12, T, label="", color="blue")
scatter!(1:10, b, label="Dados")
for i=1:10
    plot!([i, i], [b[i], T(i)], label="", color="orange", linestyle=:dashdot)
end

png(p, "plot7")