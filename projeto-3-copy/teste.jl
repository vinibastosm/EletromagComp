using Plots

# Definições iniciais
Nx = 101
Ny = 101
c = 1
dx = 1
dy = 1
dt = (min(dx, dy) / (sqrt(2))*c)
sigma = 0
epsolon_0 = 1
#epsolon_r = 1
courant = c*dt/dx

B_z = zeros(Nx, Ny)
E_x = zeros(Nx, Ny + 1)
E_y = zeros(Nx + 1, Ny)

epsolon_r = ones(Nx, Ny)

epsolon_r[:,51:end] .= 100

s = (dt*sigma)./(epsolon_r .* epsolon_0)

# Fields at time n, n-1 for Mur ABC
B_z_n = copy(B_z)
B_z_n_1 = copy(B_z_n)
E_x_n = copy(E_x)
E_x_n_1 = copy(E_x_n)
E_y_n = copy(E_y)
E_y_n_1 = copy(E_y_n)

# Mur ABC coefficients
dtdx = sqrt(dt / dx * dt / dy)
dtdx_2 = 1 / dtdx + 2 + dtdx
c_0 = -(1 / dtdx - 2 + dtdx) / dtdx_2
c_1 = -2 * (dtdx - 1 / dtdx) / dtdx_2
c_2 = 4 * (dtdx + 1 / dtdx) / dtdx_2

source_x = div(Nx,2)
source_y = div(Ny,2)

n_iter = 300

# Nome do diretório
dir_name = "imagens_heatmap"

# Verifica se o diretório existe; se não, cria
if !isdir(dir_name)
    mkdir(dir_name)
end


for n in 1:n_iter
    # Atualiza o campo magnético na etapa de tempo n+1/2
    diff_E_x = dt / dy * (E_x[:, 2:end] - E_x[:, 1:end-1])
    diff_E_y = dt / dx * (E_y[2:end, :] - E_y[1:end-1, :])
    global B_z = B_z - (courant*diff_E_y -courant*diff_E_x)

    # Loop para atualizar E_x
    for i in 2:size(E_x, 2)-1  # Percorre as colunas de E_x
        for j in 1:size(E_x, 1)  # Percorre as linhas de E_x
            E_x[j, i] = - (s[j, i] + 1) / (s[j, i] - 1) * E_x[j, i] + (courant / epsolon_r[j, i]) * (1 / (s[j, i] + 1)) * (B_z[j, i] - B_z[j, i - 1])
        end
    end

    # Loop para atualizar E_y
    for i in 2:size(E_y, 1)-1  # Percorre as linhas de E_y
        for j in 1:size(E_y, 2)  # Percorre as colunas de E_y
            E_y[i, j] = - (s[j, i] + 1) / (s[j, i] - 1) * E_y[i, j] - (courant / epsolon_r[j, i]) * (1 / (s[j, i] + 1)) * (B_z[i, j] - B_z[i - 1, j])
        end
    end

    # Pulso na etapa de tempo n+1
    tp = 30
    if n * dt <= tp
        C = (7 / 3) ^ 3 * (7 / 4) ^ 4
        pulse = C * (n * dt / tp) ^ 3 * (1 - n * dt / tp) ^ 4
    else
        pulse = 0
    end 
    
    #global B_z[source_x, source_y] = B_z[source_x, source_y] + pulse
    for i in 1:Nx
        B_z[i, 3] += pulse
    end

    # Mur ABC for left boundary
    global E_y[1, :] = c_0 * (E_y[3, :] + E_y_n_1[1, :]) + c_1 * (E_y_n[1, :] + E_y_n[3, :] - E_y[2, :] - E_y_n_1[2, :]) + c_2 * E_y_n[2, :] - E_y_n_1[3, :]

    # Mur ABC for right boundary
    global E_y[end-1, :] = c_0 * (E_y[end-3, :] + E_y_n_1[end-1, :]) + c_1 * (E_y_n[end-1, :] + E_y_n[end-3, :] - E_y[end-2, :] - E_y_n_1[end-2, :]) + c_2 * E_y_n[end-2, :] - E_y_n_1[end-3, :]

    # Mur ABC for bottom boundary
    global E_x[:, 1] = c_0 * (E_x[:, 3] + E_x_n_1[:, 1]) + c_1 * (E_x_n[:, 1] + E_x_n[:, 3] - E_x[:, 2] - E_x_n_1[:, 2]) + c_2 * E_x_n[:, 2] - E_x_n_1[:, 3]

    # Mur ABC for right boundary
    global E_x[:, end-1] = c_0 * (E_x[:, end-3] + E_x_n_1[:, end-1]) + c_1 * (E_x_n[:, end-1] + E_x_n[:, end-3] - E_x[:, end-2] - E_x_n_1[:, end-2]) + c_2 * E_x_n[:, end-2] - E_x_n_1[:, end-3]

    # Store magnetic and electric fields for ABC at time step n
    global E_x_n_1 = copy(E_x_n) # data for t = n-1
    global E_x_n = copy(E_x)     # data for t = n

    global E_y_n_1 = copy(E_y_n) # data for t = n-1
    global E_y_n = copy(E_y)     # data for t = n

    global B_z_n_1 = copy(B_z_n) # data for t = n-1
    global B_z_n = copy(B_z)     # data for t = n """

    # Cria o gráfico de calor
    heatmap(E_x, title="Gráfico de Calor de Ex, Iteração: $n", xlabel="X", ylabel="Y", color=:viridis)

    # Salva a imagem
    savefig("imagens_heatmap/heatmap_iter_$n.png")  # Salva como PNG
end
