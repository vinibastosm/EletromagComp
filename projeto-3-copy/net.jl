using Plots

# Definições iniciais
c = 3e8 

dx = 0.05e-6 #(microns)
dy = 0.05e-6 #(microns)

dt = (min(dx, dy) / sqrt(2)) / c

Nx = 101
Ny = 101

B_z = zeros(Nx, Ny)
E_x = zeros(Nx, Ny + 1)
E_y = zeros(Nx + 1, Ny)

sigma = 0 #condutividade - começar sem
epsilon = 8.8e-12 
epsilon_r = 1 #começar pelo vácuo

source_x = div(Nx, 2)
source_y = div(Ny, 2)

n_iter = 300 #tempo total

# Nome do diretório
dir_name = "imagens_heatmap1"

# Verifica se o diretório existe; se não, cria
if !isdir(dir_name)
    mkdir(dir_name)
end


for n in 1:n_iter
    # Atualiza o campo magnético na etapa de tempo n+1/2
    diff_E_x = dt / dy * (E_x[:, 2:end] - E_x[:, 1:end-1])
    diff_E_y = dt / dx * (E_y[2:end, :] - E_y[1:end-1, :])
    B_z .= B_z .- (diff_E_y .- diff_E_x)

    # Atualiza os campos elétricos na etapa de tempo n+1
    E_x[:, 2:end-1] .= E_x[:, 2:end-1] .+ dt / dy * (B_z[:, 2:end] .- B_z[:, 1:end-1])
    E_y[2:end-1, :] .= E_y[2:end-1, :] .- dt / dx * (B_z[2:end, :] .- B_z[1:end-1, :])

    # Pulso na etapa de tempo n+1
    tp = 30
    if n * dt <= tp
        #println("teste")
        C = (7 / 3) ^ 3 * (7 / 4) ^ 4
        pulse = C * (n * dt / tp) ^ 3 * (1 - n * dt / tp) ^ 4
    else
        pulse = 0
    end 
    
    for i in 1:Nx
        B_z[i, 1] += pulse
    end 

    # Cria o gráfico de calor
    heatmap(B_z, title="Gráfico de Calor de B_z, Iteração: $n", xlabel="X", ylabel="Y", color=:viridis)

    # Salva a imagem
    savefig("imagens_heatmap/heatmap_iter_$n.png")  # Salva como PNG
end
