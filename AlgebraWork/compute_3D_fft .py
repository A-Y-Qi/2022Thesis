from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from grid_setup import *
from fouriergrid_setup import *
from newtimestep import *
from angle_list import *
from plot_animation import *

b=5
test_matrix = random.randint(10, size=(b,b))
matrix_fft=fft.fft2(test_matrix)
row_fft = zeros((b,b), dtype=complex)
column_fft=zeros((b,b), dtype=complex)

#sample original matrix
sample_3D=[]
for i in range (0,b):
    w=random.randint(10, size=(b,b))
    sample_3D.append(w)
sample_3D=array(sample_3D)

#new temp matrix
z_fft=[]      
for i in range (0,b):
    z_fft.append([])
    for j in range (0,b):
        z_fft[i].append([])
        for k in range(0,b):
            z_fft[i][j].append(0)

    
#compute FFT
for k in range(0,b):
    test_matrix=sample_3D[k]
    for i in range(0,b):
        row_fft_i=fft.fft(test_matrix[i])
        row_fft[i]=row_fft_i
    for i in range(0,b):
        column_fft_i=fft.fft(row_fft[:,i])
        column_fft[:,i]=column_fft_i    
    z_fft[k]=column_fft
z_fft=array(z_fft)

for i in range(0,b):
    for j in range(0,b):
        z_fft_k=fft.fft(z_fft[:,i,j])
        z_fft[:,i,j]=array(z_fft_k)
    

sample_3D[:,0,0]

##test on Gaussian:

N = 64
x_ = (2*pi/N)*arange(-N/2.0,N/2.0)
y_ = (2*pi/N)*arange(-N/2.0,N/2.0)
z_ = (2*pi/N)*arange(-N/2.0,N/2.0)

x, y, z = meshgrid(x_, y_, z_, indexing='ij')


#k1_ = hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
#k2_ = hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
#k3_ = hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])

#k1, k2, k3 = meshgrid(k1_, k2_, k3_, indexing='ij')

k1_ = hstack([arange(-N/2.0,0),arange(0,N/2.0)])
k2_ = hstack([arange(-N/2.0,0),arange(0,N/2.0)])
k3_ = hstack([arange(-N/2.0,0),arange(0,N/2.0)])

k1, k2, k3 = meshgrid(k1_, k2_, k3_, indexing='ij')

K=k1*k1+k2*k2+k3*k3
gaussian=exp(-50*(x*x+y*y+z*z))
anaf_gaussian=power((sqrt(pi/50)),3)*exp(-pi*pi*1/50*K)

#Function with normalization 

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
            row_fft[i]=around(row_fft_i,decimals=10)
        for i in range(0,size):
            column_fft_i=1/sqrt(size)*fft.fftshift(fft.fft((row_fft[:,i])))
            column_fft[:,i]=around(column_fft_i,decimals=10)
        z_fft[k]=column_fft
        z_fft=array(z_fft)
    z_fft=array(z_fft)
    #z_fft=abs(array(z_fft))
    for i in range(0,size):
        for j in range(0,size):
            z_fft_k=1/sqrt(size)*fft.fftshift(fft.fft((z_fft[:,i,j])))
            z_fft[:,i,j]=array(around((z_fft_k),decimals=10))
    return z_fft


fft_result=fft_3D_norm(gaussian,N)

K1, K2=meshgrid(k1_, k2_)
#plot the figure for ABBD2
fig1 = plt.figure()
ax = axes3d.Axes3D(fig1)
u1=array(real(fft_result[:][5][:]))
ax.plot_surface(K1,K2,u1,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('k')
ax.set_ylabel('gaussian_fft')
plt.suptitle("Plot of fft calculated by program")


x_ = (2*pi/N)*arange(-N/2.0,N/2.0)
y_ = (2*pi/N)*arange(-N/2.0,N/2.0)
z_ = (2*pi/N)*arange(-N/2.0,N/2.0)

X, Y= meshgrid(x_, y_)
fig2 = plt.figure()
ax = axes3d.Axes3D(fig2)
u2=array(real(gaussian[5]))
ax.plot_surface(X,Y,u2,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('gaussian')
plt.suptitle("Plot of gaussian")

fig3 = plt.figure()
ax = axes3d.Axes3D(fig3)
u3=array(anaf_gaussian[5])
ax.plot_surface(K1,K2,u3,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('k')
ax.set_ylabel('gaussian_fft')
plt.suptitle("Plot of fft calculated analytically")

