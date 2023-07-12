#Dimensões & constantes
m=3 #Número de funções convexas que compõe f
n=5 #Dimensão da variável de decisão
k_max=100 #Número máximo de iterações
Θ=2*n #limite superior na metade do quadrado do diâmetro de C, que é a caixa [-e, e] nesse caso
Lfi=[sqrt(n), 1, sqrt(n)] #Lipschitz constant de cada fᵢ
Lf=sqrt(m*sum(Lfi.^2)) #Lipschitz constant de f 

x₀=randn(Float64, n)