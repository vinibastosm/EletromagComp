using OffsetArrays

c = 1 #velocidade da luz
temp_c = 1 #tempo caracteristico
Lx = 10/(temp_c*c) #comprimeNto do domínio em x  
Ly = 10/(temp_c*c) #comprimeNto do domínio em y
Nx = 10 # número de poNtos na grade em x
Ny = 10 # número de poNtos na grade em y
dx = Lx / Nx # passo do espaço grade em x
dy = Ly / Ny # passo na espaço em y
dt = 10 # passo de tempo
T = 100 # tempo total de simulção 
Nt = div(T,dt) # número total de passos de tempo
B_0 = 10 #iNtencidade campo B em T
E_0 = 10 #iNtencidade campo E em V/m

Ex =  OffsetArray(zeros(Nt+1, Nx, Ny), -1:Nt-1 ,0:Nx-1, 0:Ny-1)
Ey =  OffsetArray(zeros(Nt+1, Nx, Ny), -1:Nt-1 ,0:Nx-1, 0:Ny-1)
cBz = OffsetArray(zeros(Nt+1, Nx, Ny), -1:Nt-1 ,0:Nx-1, 0:Ny-1)

#condições iniciais para q=-1 e q=0
q_aux=0
n_aux=0
for m in 0:2:Nx-1
    println(q_aux-1,m,n_aux+1)
    Ex[q_aux-1,m,n_aux+1] = E_0
    cBz[q_aux,m,n_aux] = c*B_0
end

sy = dy*c/dt
sx = dx*c/dt

for q in 0:2:Nt-3
    for m in 0:2:Nx-3
        for n in 0:2:Ny-3
            Ex[q+1,m,n+1] = -Ex[q-1,m,n+1] + (cBz[q,m,n+2] - cBz[q,m,n])
            Ey[q+1,m+1,n] = -Ey[q-1,m+1,n] + (cBz[q,m,n] - cBz[q,m+2,n])
            if n==0 
                Ex[q+1,m,n+1] = E_0
            end

            cBz[q+2,m,n] = cBz[q,m,n] + sy*(Ex[q+1,m,n+1]-Ex[q+1,m,n-1]) - sx*(Ex[q+1,m+1,n]-Ex[q+1,m-1,n])
            if n==0 
                cBz[q+2,m,n] = c*B_0
            end
        end
    end
end

