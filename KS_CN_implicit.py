#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from scipy.linalg import toeplitz
from scipy.linalg import solve
from mpl_toolkits.mplot3d import axes3d
from numpy.linalg import inv

# Set up grid and two-solution initial data:
N = 128
dt = 0.1
h=2*pi/N
x = h*arange(-N/2.0,N/2.0)
#Setup the spectral differential matrix
col = hstack((0, 0.5*(-1.)**(arange(1, N))/tan(arange(1, N)*h/2.)))
D = toeplitz(col, -col)
col2 = hstack(((-pi**2)/(3*h**2)-(1/6), -0.5*(-1.)**(arange(1, N))/(sin(arange(1, N)*h/2.))**2))
D2 = toeplitz(col2, -col2)
D4=matmul(D2,D2)
#initial condition for ux
u = exp(-x**2)
limit=10**(-1)


# Solve PDE by and plot results by crank-Nicolson
tmax = 10
nplt = floor((tmax/25)/dt)
nmax = int(round(tmax/dt))
udata = [u]
tdata = [0]

for n in range(1,nmax+1):
    t = n*dt
    I=eye(N)
    utemp=u
    e=1
    count=0
    #Newton's Method
    while e > limit:
        M1=2*I
        M2=dt*utemp*D
        M3=dt*diag(matmul(D,utemp))
        Mdt=M1+M2+M3+dt*D2+dt*D4
        fu1=-u*matmul(D,u)*dt
        fu2=-matmul(D2,u)*dt
        fu3=-matmul(D4,u)*dt
        futemp1=-utemp*matmul(D,utemp)*dt
        futemp2=-matmul(D2,utemp)*dt
        futemp3=-matmul(D4,utemp)*dt
        f=fu1+fu2+fu3+futemp1+futemp2+futemp3-2*utemp+2*u
        norm=(sum(abs(f)**2,axis=-1))**(0.5)
        #if norm<1:
            #print(t)
        du=matmul(inv(Mdt),f)
        utemp=utemp+du
        e=max(abs(du))
        #print(e)
        count=count+1
    u=utemp
    if ((n%int(nplt))==0):
        udata.append(list(u))
        tdata.append(t)
       
#plot the figure
fig4 = plt.figure()
X, Y = meshgrid(x, tdata)
ax = axes3d.Axes3D(fig4)


u3=array(udata)


ax.plot_surface(X,Y,u3,rstride=1, cstride=1,cmap='jet', edgecolor='none')

ax.set_xlabel('x')
ax.set_ylabel('t')
plt.show()
       
