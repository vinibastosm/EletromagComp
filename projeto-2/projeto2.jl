using LinearAlgebra;
using Plots;

gs = 30;  # Grid spacing
R = 0.5;  # Radius of loop (mm)
wr = 0.1;  # Radius of wire (mm)
p = 0.1;  # Pitch of wire, centre-to-centre (mm)
N = 100;  # Number of segments in single loop of wire
n = 1; # Number of loops of wire
mu = 1;  # Magnetic susceptibility
I = -1;  # Current
C = mu*I/(4*pi);
xmin = -2.1;
xmax = 2.1;
ymin = -2.1;
ymax = 2.1;
zmin = -1.1;
zmax =  p*n*2+1.1;
x = range(xmin, xmax, length=gs);  # Positions for xrange(zmin, zmax, length=gs) 
y = range(ymin, ymax, length=gs);  # Positions for y
z = range(zmin, zmax, length=gs);  # Positions for z
Bx = zeros(gs, gs);  # x components don't change
By = zeros(gs, gs);  # y components of field matrix
Bz = zeros(gs, gs);  # z components of field matrix
norms = zeros(gs, gs);  # matrix for norms at each point
vs = zeros(4,gs);

function main()
    theta = zeros(n*N);

    println("Bem-vindo ao programa!")

    for i in 1:n*N
        theta[i] = (i - 1) * 2 * pi / N;
    end
    
    insidez = find_field(theta, R, N, wr, p);
    print(Bz);
end

# Function to do summation over all segments of wire
function find_B(pos, theta, R, N, wr, p)

    result_cross = zeros(3);
    for k in 2:n*N

        rs = [R*cos(theta[k] - pi/N),
            R*cos(theta[k] - pi/N),
            (p*(theta[k]-pi/N))/pi];

        r = pos - rs;

        dl = [R*(cos(theta[k])-cos(theta[k-1])),
            R*(sin(theta[k]-sin(theta[k-1]))),
            p/N];
        
        if norm(r) <=1.35*wr
            inwire = [0, 0, 0];
            return inwire;
        else
            result_cross += cross(dl, r)*(C/norm(r)^3)
        end

    end

    return result_cross;
end

# Calculate the magnetic field and find norms
function find_field(theta, R, N, wr, p)

    insidez = 0;

    for j in eachindex(y)
        for k in eachindex(z)

            pos = [0, y[j], z[k]];
            Bx[j, k], By[j, k], Bz[j, k] = find_B(pos, theta, R, N, wr, p);
            norms[j, k] = norm([Bx[j, k], By[j, k], Bz[j, k]]);

            if k == length(z)/2

                vs[1, j] = Bx[j, k];
                vs[2, j] = By[j, k];
                vs[3, j] = Bz[j, k];
                vs[4, j] = abs(y[j]) - R;

                if j == length(y)/2

                    insidez = Bz[j, k]
                end

            end
        end
    end

    return insidez;
end