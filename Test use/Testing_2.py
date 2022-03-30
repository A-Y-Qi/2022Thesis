from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from matplotlib import *
from w_func import *
from Kernel_d import Kernel_d
from Kernel_o import Kernel_o
from Turning_Function import *
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_2D_IFFT import *
from compute_4D_FFT import *
from compute_4D_IFFT import *
from Kernel_combined import *
from K_w_combine import *
from Saturation import *
from S_u_combine import *
from phi_prime_int import *
from Lambda_u_combine import *
from Lambda_Function import *
from Moving_Func import *
from newtimestep import *

#comment out for testing

N=32
sigma=0.6
gamma=1
kk=1
#test on K_attract
#test3=fft_3D(K_attract,N)
position=(2*pi/N)*arange(-N/2.0,N/2.0)
X,Y=meshgrid(position, position)
phi_test=6
phi_primt_test=0
dt=0.01
#testing on the w function
#test2=fft_4D(w,N)
test_angle=0
q=1

#Part that will not change over time 
K=Kernel(N,pi/4,pi/4,"r")
print("Kernel")
w=w_func_gaussian_2(N, sigma,kk,"r")
print("W finished")
K_w=K_w_combine(N,K,w)
print("K_w finished")

##Construct initial condition
sample_2D=zeros((N, N))
sample_3D=[]
for i in range (0,N):
    sample_3D.append(sample_2D)
sample_3D=array(sample_3D) #Will use for construct u0
sample_3D_2=sample_3D
sample_3D_3=sample_3D
sample_4D=[]
for i in range (0,N):
    sample_4D.append(sample_3D)
sample_4D=array(sample_4D).astype("float64")

##Construct a testing matrix
h=2*pi/N
for i in range (0,N):
    for j in range (0,N):
        for k in range(0,N):
            sample_3D[i,j,k]=exp(pi*(cos(position[j])+cos(position[k])))
u=array(sample_3D)


##plot of u0
#Ztemp=u[16]
#Z=Ztemp.astype("float64")
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot_surface(X,Y,Z,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
#plt.suptitle("Initial function u at phi="+str(position[16]))

tmax=1
time_step=int(tmax/dt)
#Animation use will use later
animation_use=[]
for i in range (0, time_step):
    animation_use.append(sample_3D)
animation_use=array(animation_use)

#started loop for time
for i in range(0,time_step):
    test_turning=Turning_2(q,K_w, u, N)
    test_turning=array(test_turning)
    saturation_func=saturation_4D(test_turning,N)
    S_u=S_u_combine(N, saturation_func, u)
    Turning=phi_prime_int_2(N, S_u)   
    lambda_func=Lambda(N,saturation_func)
    Lambda_combine=Lambda_u_combine(N,lambda_func,u)
    moving_func=Moving(N,gamma,u)
    u_new=newtimestep(u, Turning, Lambda_combine, moving_func, dt,N)
    u=real(u_new)
    animation_use[i]=u
print("time update finished")
save('N32dt001r_2', animation_use)






#test_turning=sample_4D
#test_turning=Turning_2(q,K_w, u, N)
#test_turning=array(test_turning)
#print("turning function finish")
#saturation_func=saturation_4D(test_turning,N)
#S_u=S_u_combine(N, saturation_func, u)
#Turning=phi_prime_int_2(N, S_u)
##print(Turning)

##commented
#lambda_func=Lambda(N,saturation_func)
#Lambda_combine=Lambda_u_combine(N,lambda_func,u)
#print("Lambda whole term finished")
#moving_func=Moving(N,gamma,u)
#print("Moving finished")





#plot_use3=real(Lambda_combine[31])
#fig3 = plt.figure()
#ax3 = fig3.gca(projection='3d')
#ax3.plot_surface(X,Y,plot_use3,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
#ax3.set_xlabel('x')
#ax3.set_ylabel('y')
##ax1.set_zlabel('Turning function with u=exp(pi*(cos(x)+cos(y)))')
#plt.suptitle("Lambda u function phi'="+str(position[31]))

##plot_use4=Lambda_combine[30]
##fig4 = plt.figure()
##ax4 = fig4.gca(projection='3d')
##ax4.plot_surface(X,Y,plot_use4,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
##ax4.set_xlabel('x')
##ax4.set_ylabel('y')
###ax1.set_zlabel('Turning function with u=exp(pi*(cos(x)+cos(y)))')
##plt.suptitle("Turning function with u at phi="+str(position[1])+" phi'="+str(position[0]))

##plot_use5=moving_func[0]
##fig5 = plt.figure()
##ax5 = fig5.gca(projection='3d')
##ax5.plot_surface(X,Y,plot_use5,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
##ax5.set_xlabel('x')
##ax5.set_ylabel('y')
###ax5.set_zlabel('Turning function with u=exp(pi*(cos(x)+cos(y)))')
##plt.suptitle("Turning function with u at phi="+str(position[1])+" phi'="+str(position[0]))

#plot_use6=u[7]
#fig6 = plt.figure()
#ax6 = fig6.gca(projection='3d')
#ax6.plot_surface(X,Y,plot_use6,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
#ax6.set_xlabel('x')
#ax6.set_ylabel('y')
##ax1.set_zlabel('Turning function with u=exp(pi*(cos(x)+cos(y)))')
#plt.suptitle("Turning function with u at phi="+str(position[1])+" phi'="+str(position[0]))





