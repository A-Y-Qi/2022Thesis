from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy.linalg import *
from compute_2D_FFT import *
from compute_2D_IFFT import *
from Matrix_A import *
from Matrix_B import *
from matplotlib import *
from matplotlib import animation
from Arrow_function import *
#from Arrow_function import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (k and l)
## where [I-dtM]^(-1)

N=64
position=(2*pi/N)*arange(-N/2.0,N/2.0)
X,Y=meshgrid(position, position)
dt=0.01
u=load("u_file_name.npy")
sum_phi=load("sum_phi_file_name.npy")
Ex_t,Ey_t, x_t,y_t=Arrow_max1000(u,sum_phi,N)
timestep=len(u)

            
'''
Plot sum phi
'''
time=0  
fig1 = plt.figure()
ax1 = fig1.gca()
#ax1 = fig.gca()
Z1 = sum_phi[time]
con1=ax1.contourf(X, Y, Z1,cmap=plt.cm.terrain)
ax1.quiver(x_t[time],y_t[time],Ex_t[time],Ey_t[time], edgecolors='w',color='w')
fig1.colorbar(con1)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
plt.suptitle("the Gaussian function with saturation term at t="+str(time*dt))

time=50
fig2 = plt.figure()
ax2 = fig2.gca()
Z2 = sum_phi[time]
con2=ax2.contourf(X, Y, Z2,cmap=plt.cm.terrain)
ax2.quiver(x_t[time],y_t[time],Ex_t[time],Ey_t[time], edgecolors='w',color='w')
fig2.colorbar(con2)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
plt.suptitle("the Gaussian function with saturation term at t="+str(time*dt))

time=100
fig3 = plt.figure()
ax3 = fig3.gca()
Z3 = sum_phi[time]
con3=ax3.contourf(X, Y, Z3,cmap=plt.cm.terrain)
ax3.quiver(x_t[time],y_t[time],Ex_t[time],Ey_t[time], edgecolors='w',color='w')
fig3.colorbar(con3)
ax3.set_xlabel('x')
ax3.set_ylabel('y')
##ax6.view_init(90, 90)
plt.suptitle("the Gaussian function with saturation term at t="+str(time*dt))

time=150
fig4 = plt.figure()
ax4 = fig4.gca()
Z4 = sum_phi[time]
con4=ax4.contourf(X, Y, Z4,cmap=plt.cm.terrain)
ax4.quiver(x_t[time],y_t[time],Ex_t[time],Ey_t[time], edgecolors='w',color='w')
fig4.colorbar(con4)
ax4.set_xlabel('x')
ax4.set_ylabel('y')
plt.suptitle("the Gaussian function with saturation term at t="+str(time*dt))

time=200
fig5 = plt.figure()
ax5 = fig5.gca()
Z5 = sum_phi_r[time]
con5=ax5.contourf(X, Y, Z5,cmap=plt.cm.terrain)
ax5.quiver(x_t[time],y_t[time],Ex_t[time],Ey_t[time], edgecolors='w',color='w')
fig5.colorbar(con5)
ax5.set_xlabel('x')
ax5.set_ylabel('y')
plt.suptitle("the Gaussian function with saturation term at t="+str(time*dt))

time=250
fig6 = plt.figure()
ax6 = fig6.gca()
Z6 = sum_phi_r[time]
con6=ax6.contourf(X, Y, Z6,cmap=plt.cm.terrain)
ax6.quiver(x_t[time],y_t[time],Ex_t[time],Ey_t[time], edgecolors='w',color='w')
fig6.colorbar(con6)
ax6.set_xlabel('x')
ax6.set_ylabel('y')
plt.suptitle("the Gaussian function with saturation term at t="+str(time*dt))


'''
Animation
'''
fig7 = plt.figure()
ax7 = fig7.gca(projection='3d')
#ax1 = fig.gca()
#timestep_number=99

kk = 0
def animate2(i):
    global kk
    Z7 = sum_phi[kk]
    kk += 1
    ax7.clear()
    ax7.plot_surface(X,Y,Z7,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
    ax7.set_zlim(0,1000)
    ax7.set_xlabel('x')
    ax7.set_ylabel('y')
    #ax1.view_init(90, 90)
    #plt.suptitle("Simulation using Matrix when N="+str(N)+" phi="+str(position[kk])+" dt="+str(dt))
    plt.suptitle("Simulation of ...... at t="+str(kk*dt)[:4])
    if kk >= timestep:
        kk=0

ax7.view_init(90, 90)
anim7 = animation.FuncAnimation(fig7,animate2)
plt.show()  
