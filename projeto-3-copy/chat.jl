using PyPlot
using Printf
using ImageMagick
using FileIO

#inciando o problema
c = 1 

dx = 1 #(microns)
dy = 1 #(microns)

dt = dx/(c*sqrt(2)) #(ordem de femtosegundos) e courant = 1/sqrt()

println(dt)
println((c*dt)/dx)

Ns = 50
Nt = 100
Tempo = dt * Nt
Comp = dx * Ns 

sigma = 10e-10 #condutividade - começar sem
epsilon = 8.8e-12 
epsilon_r = 1 #começar pelo vácuo

courant = 1/sqrt(2) # = (c*dt)/dx
println(courant)
s = (dt*sigma)/(epsilon*epsilon_r)

# inicializando os campos com zeros
# (x, y, tempo)
Ex = zeros(Ns+1, Ns+1, Nt+1)
Ey = zeros(Ns+1, Ns+1, Nt+1)
cBz = zeros(Ns+1, Ns+1, Nt+1)


#loop temporal

# se q for par + 1 (impar) -> calculo cBz
# cBz(m,n,q+2) = cBz(m,n,q) + courant*(Ex(m,n+1,q+1) - Ex(n,n-1,q+1)) - courant*(Ey(m+1,n,q+1) - Ey(m-1,n,q+1)) 

# se q for impar + 1  (par) -> calculo E

# se m par + 1 (impar) e n impar + 1 (par) -> Ex
# Ex(m,n+1,q+1) = (courant/epsilon_r)*(1/(s+1))*(cBz(m,n+2,q) - cBz(m,n,q)) - ((s-1)/(s+1))*Ex(m,n+1,q-1)

# se m impar + 1 (par) e n par + 1 (impar) -> Ex
# Ey(m+1,n,q+1) = (courant/(epsilon_r*(s+1)))*(cBz(m+2,n,q) - cBz(m,n,q)) - ((s-1)/(s+1))*Ey(m+1,n,q-1)

cBz[div(Ns,2)+1,div(Ns,2)+1,3] = 10

for q in 1:Nt-1

    if q % 2 == 0
        #println(q)
        for m in 1:Ns-1
            for n in 1:Ns-1
                if m % 2 == 0 && n % 2 == 1
                    #println(m,n)
                    global Ey[m+1,n,q+1] = (courant/(epsilon_r*(s+1)))*(cBz[m+2,n,q] - cBz[m,n,q]) - ((s-1)/(s+1))*Ey[m+1,n,q-1]
                

                elseif m % 2 == 1 && n % 2 == 0
                    #println(m,n)
                    global Ex[m,n+1,q+1] = (courant/epsilon_r)*(1/(s+1))*(cBz[m,n+2,q] - cBz[m,n,q]) - ((s-1)/(s+1))*Ex[m,n+1,q-1]

                end
            end
        end

    
    elseif q % 2 == 1
        #println(q)
        for m in 3:2:Ns-1
            for n in 3:2:Ns-1
                #println(m,n)
                global cBz[m,n,q+2] = cBz[m,n,q] + courant*(Ex[m,n+1,q+1] - Ex[m,n-1,q+1]) - courant*(Ey[m+1,n,q+1] - Ey[m-1,n,q+1])
            end
        end
       
    end
end
# Função para plotar o campo Bz no tempo q e salvar a imagem
function plot_bz_field(cBz, q)
    figure()
    imshow(cBz[:, :, q]', cmap="inferno", origin="lower", extent=(0, Ns*dx*1e6, 0, Ns*dy*1e6))
    colorbar(label="Bz")
    title("Campo Magnético Bz no tempo q = $q")
    xlabel("Posição x (µm)")
    ylabel("Posição y (µm)")
    filename = @sprintf("bz_field_q%03d.png", q)
    savefig(filename)
    close()
    return filename
end

# Gerando as imagens e armazenando os nomes dos arquivos
frames = []
for q in 1:Nt
    push!(frames, plot_bz_field(cBz, q))
end

# Carregar as imagens e criar o GIF
images = [load(frame) for frame in frames]
gif(images, "campo_magnetico_evolucao_pulso.gif", fps=10)
println("GIF salvo em: campo_magnetico_evolucao_pulso.gif")
