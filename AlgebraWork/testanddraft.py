from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from plot_animation import *
from scipy.linalg import dft
    
    

def fft_3D_norm(matrixin3D,size):
    row_fft = zeros((size,size), dtype=complex)
    column_fft=zeros((size,size), dtype=complex)
    z_fft=[]      
    for i in range (0,size):
        z_fft.append([])
        for j in range (0,size):
            z_fft[i].append([])
            for k in range(0,size):
                z_fft[i][j].append(1) 
    #fft loop 
    for k in range(0,size):
        test_matrix=matrixin3D[k]
        for i in range(0,size):
            row_fft_i=1/sqrt(size)*fft.fftshift(fft.fft((test_matrix[i])))
            row_fft[i]=(row_fft_i)
        for i in range(0,size):
            column_fft_i=1/sqrt(size)*fft.fftshift(fft.fft((row_fft[:,i])))
            column_fft[:,i]=(column_fft_i)
        z_fft[k]=column_fft
        z_fft=array(z_fft)
    z_fft=array(z_fft)
    #z_fft=abs(array(z_fft))
    for i in range(0,size):
        for j in range(0,size):
            z_fft_k=1/sqrt(size)*fft.fftshift(fft.fft((z_fft[:,i,j])))
            z_fft[:,i,j]=array((z_fft_k))
    return z_fft

N = 128
x_ = (2*pi/N)*arange(-N/2.0,N/2.0)
y_ = (2*pi/N)*arange(-N/2.0,N/2.0)
z_ = (2*pi/N)*arange(-N/2.0,N/2.0)

x, y, z = meshgrid(x_, y_, z_, indexing='ij')
gaussian=exp(-50*x*x-y*y-z*z)





fft_result=fft_3D_norm(gaussian,N)


#test on normalizaion of analytic solution
k1_ = hstack([arange(-N/2.0,0),arange(0,N/2.0)])
k2_ = hstack([arange(-N/2.0,0),arange(0,N/2.0)])
k3_ = hstack([arange(-N/2.0,0),arange(0,N/2.0)])

k1, k2, k3 = meshgrid(k1_, k2_, k3_, indexing='ij')

K=k1*k1+k2*k2+k3*k3

anaf_gaussian1=power((sqrt(pi/50)),3)*exp(-1/50*pi*pi*K)

plot_gaussianF_1=anaf_gaussian1[5]

K1,K2=meshgrid(k1_, k2_)
fig1 = plt.figure()
ax = axes3d.Axes3D(fig1)
ax.plot_surface(K1,K2,plot_gaussianF_1,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('k')
ax.set_ylabel('gaussian_fft')
plt.suptitle("Plot of fft calculated by continuous fourier")





def dft_3D_norm(matrixin3D,size):
    m = dft(size)
    row_fft = zeros((size,size), dtype=complex)
    column_fft=zeros((size,size), dtype=complex)
    z_fft=[]      
    for i in range (0,size):
        z_fft.append([])
        for j in range (0,size):
            z_fft[i].append([])
            for k in range(0,size):
                z_fft[i][j].append(1) 
    #Dft loop 
    for k in range(0,size):
        test_matrix=matrixin3D[k]
        for i in range(0,size):
            row_fft_i=1/sqrt(size)*fft.fftshift(m @ test_matrix[i])
            row_fft[i]=row_fft_i
        for i in range(0,size):
            column_fft_i=1/sqrt(size)*fft.fftshift(m @ row_fft[:,i])
            column_fft[:,i]=(column_fft_i)
        z_fft[k]=column_fft
        z_fft=array(z_fft)
    z_fft=array(z_fft)
    #z_fft=abs(array(z_fft))
    for i in range(0,size):
        for j in range(0,size):
            z_fft_k=1/sqrt(size)*fft.fftshift(m @ z_fft[:,i,j])
            z_fft[:,i,j]=array(abs(z_fft_k))
    return z_fft
        
gaussianF_DFT=dft_3D_norm(gaussian,N)
plot_gaussianF_DFT=real(gaussianF_DFT[:][5][:])
fig2 = plt.figure()
ax = axes3d.Axes3D(fig2)
K11,K22=meshgrid(k1_, k2_)
ax.plot_surface(K11,K22,plot_gaussianF_DFT,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('k')
ax.set_ylabel('gaussian_fft')
plt.suptitle("Plot of Gaussian Fourier calculated by DFT")

M=dft(N)
A=fft.fft((gaussian[5][2][:]))
B=M @ gaussian[5][2][:]
c=linalg.norm(A-B)