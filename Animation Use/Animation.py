from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from matplotlib import *

N=32
position=(2*pi/N)*arange(-N/2.0,N/2.0)
X,Y=meshgrid(position, position)
dt=0.01
animation_use=load("N32dt001r_2.npy")
###Animation use
fig = plt.figure()
ax1 = fig.gca(projection='3d')
#ax1 = fig.gca()
timestep_number=99

k = 0
def animate(i):
    global k
    ##k is time and 9 is the angel position [9]
    Z = animation_use[k,9]
    k += 1
    ax1.clear()
    ax1.plot_surface(X,Y,Z,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
    #ax1.contourf(X,Y,Z)
    ax1.set_zlim(0,400)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    #ax1.view_init(90, 90)
    plt.suptitle("Simulation using Euler when N="+str(N)+" t="+str(k*0.01)+" dt="+str(dt))
    if k > timestep_number-1:
        k=0

    
anim = animation.FuncAnimation(fig,animate)
plt.show()    

#animation_use_1=load("N32dt001r.npy")
#animation_use_2=load("N32dt001a.npy")
#diff=animation_use_1-animation_use_2
#max_num=0
#ii=0
#jj=0
#kk=0
#for i in range (0,99):
    #for j in range (0,N):
        #for k in range (0,N):
            #max_temp=max(abs(diff[i,j,k]))
            #if max_temp>max_num:
                #max_num=max_temp
                #ii=i
                #jj=j
                #kk=k
#print(max_num)
