from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from w_func import *
from Kernel_d import Kernel_d
from Kernel_o import Kernel_o
from Turning_Function import *
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_4D_FFT import *
from compute_4D_IFFT import *
from Kernel_combined import *
from K_w_combine import *

#comment out for testing

N=32
sigma=pi/4
kk=0.5
K=Kernel(N,pi/4,pi/4,"a")
print("Kernel")
#test on K_attract
#test3=fft_3D(K_attract,N)
position=(2*pi/N)*arange(-N/2.0,N/2.0)
X,Y=meshgrid(position, position)


w=w_func_gaussian(N, sigma,kk,"a")
print("W finished")
#testing on the w function
#test2=fft_4D(w,N)


q=1
sample_2D=zeros((N, N),dtype=complex)
sample_3D=[]
for i in range (0,N):
    sample_3D.append(sample_2D)
sample_3D=array(sample_3D)
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

#plot of u0
#Ztemp=u[1]
#Z=Ztemp.astype("float64")
#fig = plt.figure()
#ax1 = fig.gca(projection='3d')
#ax1.plot_surface(X,Y,Z,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)


K_w=K_w_combine(N,K,w)
print("K_w printed")
test_turning=sample_4D
test_turning=Turning_2(q,K_w, u, N)
test_turning=array(test_turning)
#for i in range (0,N):
    #for j in range (0,N):
        #for k in range(0,N):
            #for l in range (0,N):
                #if test_FFT_turning[i,j,k,l] >1:
                    #print ("the value of fft is"+str(test_FFT_turning[i,j,k,l]) + " at "+"i,j,k,l= "+str(i)+","+str(j)+","+str(k)+","+str(l)+".")

#plot of fft_u
#test_u=fft_3D(u,N)
#test_u_2=abs(real(test_u[int(N/2)]))
#Z2temp=test_u_2
#Z2=Z2temp.astype("float64")
#fig2 = plt.figure()
#ax2 = fig2.gca(projection='3d')
#ax2.plot_surface(X,Y,Z2,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)


#test use
#max_n=0
#for i in range (0,N):
    #for j in range (0,N):
        #for k in range(0,N):
            #temp_max=max(abs(turning[i,j,k]))
            #if temp_max>max_n:
                #max_n=temp_max
#print(max_n)
plot_use=real(test_turning[0,9])
fig = plt.figure()
ax1 = fig.gca(projection='3d')
ax1.plot_surface(X,Y,plot_use,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
#ax1.set_zlabel('Turning function with u=exp(pi*(cos(x)+cos(y)))')
plt.suptitle("Turning function at phi="+str(position[0])+" phi'="+str(position[0]))






