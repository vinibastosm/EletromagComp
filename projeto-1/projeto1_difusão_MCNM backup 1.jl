using LinearAlgebra
using Plots
using LaTeXStrings

# Parâmetros
L = 10.0  # comprimento da barra
nx = 100  # número de pontos espaciais
dx = L / (nx - 1)  # passo espacial
alpha = 1.0  # difusividade térmica
k = 1.0  # condutividade térmica
a = 1.0  # raio da fonte
P0 = 10.0  # potência máxima
ω = 2.0 * π  # frequência angular
xf = 5.0  # posição da fonte
nt = 500  # número de passos de tempo
dt = 0.01  # passo de tempo
T_inf = 1.0  # temperatura do infinito
T_0 = 1.0  # temperatura fixa da fonte

# Identificar o índice da posição da fonte na malha espacial
source_index = Int(floor(xf / dx))

# Criação da malha espacial e inicialização da temperatura
x = range(0, stop=L, length=nx)
T = fill(T_inf, nx)  # Temperatura inicial constante (T_inf)
T_new = zeros(nx)

# Parâmetros do método de Crank-Nicolson
r = alpha * dt / (2 * dx^2)

# Matrizes tridiagonais
A = zeros(nx, nx)
B = zeros(nx, nx)

# Montagem das matrizes A e B
for i in 2:nx-1
    A[i, i-1] = -r
    A[i, i] = 1 + 2*r
    A[i, i+1] = -r

    B[i, i-1] = r
    B[i, i] = 1 - 2*r
    B[i, i+1] = r
end

# Condições de contorno: fronteiras esquerda e direita isoladas
A[1, 1] = 1
B[1, 1] = 1
A[end, end] = 1
B[end, end] = 1

# Função que calcula o termo fonte
P(t) = P0 * sin(ω * t)^2

# Função para atualizar a temperatura na posição da fonte
function aplicar_fluxo_fonte!(T, t)
    dT = P(t) / (k * 4 * π * a^2) * dx
    T[source_index - 1] = T_0 + dT
    T[source_index + 1] = T_0 + dT
end

# Listas para armazenar as soluções e os pontos de maior temperatura ao longo do tempo
T_list = []
x_max_left = []
x_max_right = []
v_max_left = []
v_max_right = []
last_x_max_left = 0.0
last_x_max_right = 0.0

# Solução no tempo
for n in 1:nt
    # Montar o vetor do lado direito
    b = B * T

    # Aplicar condição da fonte de calor
    aplicar_fluxo_fonte!(b, n * dt)

    # Resolução do sistema linear A * T_new = b
    global T_new = A \ b  # Declare `T_new` como global

    # Garantir que a temperatura da fonte permaneça fixa em T_0
    T_new[source_index] = T_0

    # Atualizar a solução
    T .= T_new

    # Armazenar a solução a cada passo de tempo
    push!(T_list, copy(T))

    # Encontrar os pontos de maior temperatura à esquerda e à direita da fonte
    max_left_index = argmax(T[1:source_index])
    max_right_index = argmax(T[source_index+1:end]) + source_index

    push!(x_max_left, max_left_index * dx)
    push!(x_max_right, max_right_index * dx)

    # Calcular a velocidade dos pontos de máxima temperatura
    if n > 1
        global v_max_left, v_max_right, last_x_max_left, last_x_max_right  # Declare todas as variáveis globais

        push!(v_max_left, (x_max_left[end] - last_x_max_left) / dt)
        push!(v_max_right, (x_max_right[end] - last_x_max_right) / dt)
    else
        push!(v_max_left, 0.0)
        push!(v_max_right, 0.0)
    end

    last_x_max_left = x_max_left[end]
    last_x_max_right = x_max_right[end]
end


# Função para atualizar o gráfico
function update_plot(n)
    plot(x, T_list[n], label="Temperatura")

    # Adicionar linha para T_0 (temperatura da fonte)
    hline!([T_0], color=:red, linestyle=:dash, label="T_0 = $T_0")

    # Adicionar linha para T_inf (temperatura do infinito)
    hline!([T_inf], color=:blue, linestyle=:dash, label="T_inf = $T_inf")

    # Marcar os pontos de maior temperatura à esquerda e à direita da fonte
    scatter!([x_max_left[n]], [T_list[n][Int(round(x_max_left[n] / dx))]], color=:green, label="T_max Left")
    scatter!([x_max_right[n]], [T_list[n][Int(round(x_max_right[n] / dx))]], color=:magenta, label="T_max Right")

    xlabel!("x (m)")
    ylabel!("Temperatura (T)")
    title!("Tempo = $(n*dt) s")
    
    # Definir os limites do eixo y
    ylims!(T_inf - 0.1, maximum(T_list[n]) + 0.1)
end

# Criar a animação
anim = @animate for n in 1:nt-1
    update_plot(n)
end

# Salvar animação
gif(anim, "heat_diffusion.gif", fps=20)

# Ajustar o tempo para coincidir com o número de elementos em x_max_left e x_max_right
time_positions = (1:length(x_max_left)) .* dt  # Corrigindo o vetor de tempo

# Posições dos pontos de maior temperatura
plot(time_positions, x_max_left, label="x_max Left", color=:green)
plot!(time_positions, x_max_right, label="x_max Right", color=:magenta, xlabel="Tempo (s)", ylabel="Posição x_max (m)")

# Velocidades dos pontos de maior temperatura
plot(time_positions, v_max_left, label="v_max Left", color=:green)
plot!(time_positions, v_max_right, label="v_max Right", color=:magenta, xlabel="Tempo (s)", ylabel="Velocidade v_max (m/s)")
