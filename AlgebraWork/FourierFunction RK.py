#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from grid_setup import *
from fouriergrid_setup import *
from newtimestep import *
from angle_list import *
from plot_animation import *

# Set up initial data:
N = 64
m = 64
dt = 0.01
g1=0.2 
g2=0.9
G=g1+0.5*g2
gamma=1

#Set up spacial grid
[xx,yy,X,Y]=grid_setup(N)

#initial consition with constant angle
u1=[exp(cos(pi*yy))]*m
#u_hat=fft.fft2(u1) #2D fourier transform 
###Do I need to shift in this case? 
##Alt thought 
u1_hat=fft.fft(exp(cos(pi*yy)))
u_hat=[u1_hat]*m

#Probability function 
mu, sigma = 0, 0.5
w=random.normal(mu, sigma, m)
w=w*(m/(sum(w)*2*pi))
w_hat=real(abs(fft.fft(random.normal(mu, sigma, m))))/m
F_list=-G+G*2*pi/N*w_hat


#set up fourier grid
[kk, ll, plus, minus]=fouriergrid_setup(N)


#used for plot
udata1 = [u1]
unew_hat=u_hat

##Loop starts
tmax = 15
nmax = int(round(tmax/dt))

for t in range(0,nmax):
    for q in range(0, m):
        F=F_list[q]
        p=q-1
        r=q+1
        if q==0:
            p=m-1
        if q==m-1:
            r=1
        uq_hat=u_hat[q]
        ur_hat=u_hat[r]
        up_hat=u_hat[p]
        unew_hat[q]=newtimestepRK2(uq_hat, ur_hat, up_hat, F, gamma, dt, plus, minus)
        u1[q] = real(fft.ifft(unew_hat[q]))
    #make a 3D array at different time 
    udata1.append(list(u1))
    u_hat=unew_hat

##plot the density at specific angle
    #plot_animation(udata1,0,N,X,Y)


##plot test

#fig1 = plt.figure()
#ax = axes3d.Axes3D(fig1)
#u1=array(angle_list(0,udata1,N))
#ax.plot_surface(X,Y,u1[17],rstride=1, cstride=1,cmap='viridis', edgecolor='none')
#ax.set_xlabel('x')
#ax.set_ylabel('t')
#plt.suptitle("Simulation of convolution using AB/BDI2 with sigma="+str(sigma)+" dt="+str(dt)+"N="+str(N))    


#Animation test

density=array(angle_list(0,udata1,N))

fig = plt.figure()
ax1 = fig.gca(projection='3d')
#ax1 = fig.gca()

k = 0
def animate(i):
    global k
    Z = density[k]
    k += 30
    ax1.clear()
    ax1.plot_surface(X,Y,Z,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
    #ax1.contourf(X,Y,Z)
    ax1.set_zlim(0,2.5)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.suptitle("Simulation using RK2 when N="+str(N)+" m="+str(m)+" dt="+str(dt))
    if k > nmax-1:
        k=0

    
anim = animation.FuncAnimation(fig,animate)
#ax1.view_init(90, 90)
plt.show()    
    



        


