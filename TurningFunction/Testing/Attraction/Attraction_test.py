from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from w_func import w_func
from Kernel_d import Kernel_d
from Kernel_o import Kernel_o
from Fourier_Turning_Function import Fourier_Turning
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_4D_FFT import *
from compute_4D_IFFT import *
from Kernel_attraction import *

#comment out for testing

N=32
sigma=pi/2
k=0.5
K_attract=Kernel_attraction(N,pi/4,pi/4)
#test on K_attract
#test3=fft_3D(K_attract,N)
#position=(2*pi/N)*arange(-N/2.0,N/2.0)
#X,Y=meshgrid(position, position)


w=w_func(N, sigma,k)
#testing on the w function
#test2=fft_4D(w,N)


q=0.5
sample_2D=zeros((N, N),dtype=complex)
sample_3D=[]
for i in range (0,N):
    sample_3D.append(sample_2D)
sample_3D=array(sample_3D)
sample_3D_2=sample_3D
sample_3D_3=sample_3D

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



test_FFT_turning=Fourier_Turning(0.5,K_attract,w, u, N)
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

Turning_comp=ifft_4D(test_FFT_turning, N)
turning=Turning_comp.astype("float64")
print(turning)





