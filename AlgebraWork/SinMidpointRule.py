from numpy import *
from MidpointRule import *
import matplotlib.pyplot as plt
from scipy.linalg import dft


N = 256
k = (2*pi/N)*arange(-N/2.0,N/2.0) #wave number 

fourier = lambda y: 2*cos(y)*exp(-2*pi*1j*k*y)

F_hat=midpoint(fourier, -pi, pi, N )



function=cos(k)


##for DFT
result=zeros(N, dtype=complex)
#for i in range(0,N-1):
    #a=0
    #for j in range(0,N-1):
        #a=a+function[j]*exp(-2*pi*1j*i*j/N)
    #result[i]=a

#fast=fft.fft(function)

##From paper:
lamb=1/2
for s in range(0,N-1):
    add=0
    for p in range (-N,N):
        nv_2=1-s+p*N
        if nv_2==0 or nv_2==2:
            print(nv_2)
            add=add+exp(-pi*1j*(nv_2-1))
    result[s]=add
#still troublesome
        

    
    
    


difference=F_hat-result
#print(difference)

plt.plot(k, real(result))
plt.plot(k, real(F_hat))

