import scipy.constants as constants
import matplotlib.pyplot as plt
from scipy import optimize
import matplotlib as mpl
import numpy as np
from matplotlib.colors import SymLogNorm
from PIL import Image
from io import BytesIO  # Para salvar as imagens em memória

# Constantes
eps = constants.epsilon_0
pi = constants.pi
e = constants.e
c = constants.c
u_0 = constants.mu_0

# Caracteristicas da carga 
pos_charge= True
initial_speed = 0.6*c
x_stop = 30e-9
deceleration = 0.5*initial_speed**2/x_stop
stop_t = initial_speed/deceleration

def xpos(t):
    return initial_speed*t

def ypos(t):
    return 0

def zpos(t):
    return 0

def xvel(t):
    return initial_speed

def yvel(t):
    return 0

def zvel(t):
    return 0

def xacc(t):
    return 0

def yacc(t):
    return 0

def zacc(t):
    return 0

def retarded_time(tr, t, X, Y, Z):
    return ((X-xpos(tr))**2 + (Y-ypos(tr))**2 + (Z-zpos(tr))**2)**0.5 - c*(t-tr)

def _calculate_individual_E(tr, X, Y, Z):
    rx = X - xpos(tr)
    ry = Y - ypos(tr)
    rz = Z - zpos(tr)
    r_mag = (rx**2 + ry**2 + rz**2)**0.5
    vx = xvel(tr) 
    vy = yvel(tr)
    vz = zvel(tr)
    ax = xacc(tr)  
    ay = yacc(tr)
    az = zacc(tr)
    ux = c*rx/r_mag - vx  
    uy = c*ry/r_mag - vy
    uz = c*rz/r_mag - vz
    r_dot_u = rx*ux + ry*uy + rz*uz
    r_dot_a = rx*ax + ry*ay + rz*az
    vel_mag = (vx**2 + vy**2 + vz**2)**0.5
    const = e/(4*pi*eps) * r_mag/(r_dot_u)**3
    if not pos_charge:  
        const *= -1
    xvel_field = const*(c**2-vel_mag**2)*ux
    yvel_field = const*(c**2-vel_mag**2)*uy
    zvel_field = const*(c**2-vel_mag**2)*uz
    xacc_field = const*(r_dot_a*ux - r_dot_u*ax)
    yacc_field = const*(r_dot_a*uy - r_dot_u*ay)
    zacc_field = const*(r_dot_a*uz - r_dot_u*az)

    return (xvel_field+xacc_field, yvel_field+yacc_field,zvel_field+zacc_field)

def calculate_E(t, X, Y, Z):

    t_array = np.ones((X.shape))
    t_array[:, :, :] = t
    Ex = np.zeros((X.shape))
    Ey = np.zeros((X.shape))
    Ez = np.zeros((X.shape))
    tr = optimize.newton(func=retarded_time, x0=t_array, args=(t_array, X, Y, Z), tol=1e-20)
    E_field = _calculate_individual_E(tr, X, Y, Z)
    Ex += E_field[0]
    Ey += E_field[1]
    Ez += E_field[2]

    if X.shape[0] == 1:
        return(Ex[0, :, :], Ey[0, :, :], Ez[0, :, :])
    elif X.shape[1] == 1:
        return(Ex[:, 0, :], Ey[:, 0, :], Ez[:, 0, :])
    elif X.shape[2] == 1:
        return(Ex[:, :, 0], Ey[:, :, 0], Ez[:, :, 0])

lim = 50e-9
grid_size = 1000
X, Y, Z = np.meshgrid(np.linspace(-lim, lim, grid_size), 0, np.linspace(-lim, lim, grid_size), indexing='ij')

# Calcule os valores de E_total
aux = 6
E_total = [None] * aux
pos = []
for i in range(aux):
    E_total[i] = calculate_E(t=(i) * stop_t / aux, X=X, Y=Y, Z=Z)
    pos.append(xpos((i) * stop_t / aux))

frames = []

for i in range(aux):
    # Criar o gráfico
    fig, ax = plt.subplots()
    im = ax.imshow(
        E_total[i][2].T,
        cmap='coolwarm',
        origin='lower',
        aspect='equal',
        extent=(-lim, lim, -lim, lim)
    )
    im.set_norm(SymLogNorm(linthresh=1e5, linscale=1, vmin=-1e7, vmax=1e7, base=10))
    plt.colorbar(im, label='$E_z$ [N/C]')
    plt.xlabel('$x$ [m]')
    plt.ylabel('$z$ [m]')

    # Adicionar a bolinha preta na posição correspondente
    pos_z = 0
    pos_x = pos[i]
    ax.scatter(pos_x, pos_z, color='black', s=50)

    # Salvar o gráfico em um buffer de memória
    buf = BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)

    # Converter para PIL Image enquanto o buffer ainda está aberto
    frames.append(Image.open(buf).copy())  # .copy() garante que os dados sejam mantidos
    buf.close()
    plt.close(fig)

# Criar o GIF
frames[0].save(
    "campo_eletrico.gif",
    save_all=True,
    append_images=frames[1:],
    duration=200,  # Duração de cada quadro em milissegundos
    loop=0
)

print("GIF criado com sucesso!")