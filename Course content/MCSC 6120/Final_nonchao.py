#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA

# Set up grid and two-solution initial data:
N = 256
dt = 0.05
x = (2*pi/N)*arange(-N/2.0,N/2.0)
h=2*pi/N
#initial condition for ux
u1 = exp(-x**2)
u_hat = fft.fft(u1)
uold = exp(-(x-0.2*dt)**2)
uold_hat = fft.fft(uold)
k = hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
sigma=1*10**(-1)
r=1

# Solve PDE by and plot results by only AB/BDII
tmax =20
nplt = floor((tmax/25)/dt)
nmax = int(round(tmax/dt))
udata1 = [u1]
tdata = [0]
a0=3/2*ones(len(k))
a1=-2
a2=1/2
b0=2
b1=-1
count=0

for n in range(1,nmax+1):
    t = n*dt
    g = -1j*dt*k
    #the fourier transform of u1^2
    usq_hat=fft.fft(u1**2)
    uoldsq_hat=fft.fft(uold**2)
    delta=1/sqrt(2*pi*sigma)*exp(-x**2/(2*sigma))
    G_hat=k**2*exp(-sigma*k**2/2)
    #G_hat=k**2*abs(real(fft.fft(delta)))*2*pi/N
    a=g*b0*usq_hat
    b=g*b1*uoldsq_hat
    c=a1*u_hat
    d=a2*uold_hat
    e=1/(a0-dt*(G_hat-k**4))
    unew_hat=(a+b-c-d)*e
    uold_hat=u_hat 
    u_hat=unew_hat
    u1 = real(fft.ifft(u_hat))
    uold=real(fft.ifft(uold_hat))
    count=count+1
    #if LA.norm(u1)<0.1: 
        #print(count)
    print(t,LA.norm(u1))

    if ((n%int(nplt))==0):
        udata1.append(list(u1))
        tdata.append(t)
        

#plot the figure for ABBD2
fig1 = plt.figure()
X, Y = meshgrid(x, tdata)
ax = axes3d.Axes3D(fig1)
u1=array(udata1)
ax.plot_surface(X,Y,u1,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.suptitle("Simulation with using AB/BDI2 with sigma="+str(sigma)+" dt="+str(dt)+"N="+str(N))
