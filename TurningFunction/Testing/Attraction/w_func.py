from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_4D_FFT import *
from compute_4D_IFFT import *

def w_func(N, sigma,k):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0) 
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=[]
    for i in range (0,N):
        array(w_func_3D.append(w_func_2D))  
    w_func_4D=[]
    for i in range (0,N):
        array(w_func_4D.append(w_func_3D))
    w_func_4D=array(w_func_4D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    if Sx[k]==Sy[j] and Sx[k]==0:
                        psi=0   
                    else:
                        psi=arccos(Sx[k]/(sqrt((Sx[k]**2)+(Sy[j]**2))))
                    v=k*(phi_prime[j]-psi)
                    theta=phi_prime[j]-phi[i]-v
                    if abs(theta)<=sigma:
                        w_func_4D[i,j,k,l]=1/(2*sigma)
                    else:
                        w_func_4D[i,j,k,l]=0
    return w_func_4D
