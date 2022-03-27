from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special

#Orentation kernel
#Attraction: Kernel_o_func[i,j,k]=(1/(2*pi))*(-cos(phi[i]-psi)+1)
#Repulsion: Kernel_o_func[i,j,k]=(1/(2*pi))*(cos(phi[i]-psi)+1)


def Kernel_o(N,a_or_r):
    #define theta and phi
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0)    
    #Generate a 2*2 matrix
    Kernel_o_func_1 = zeros((N,N))
    Kernel_o_func=[]
    for i in range (0,N):
        array(Kernel_o_func.append(Kernel_o_func_1))
    Kernel_o_func=array(Kernel_o_func)
    for i in range (0,N):
        for j in range(0,N):
            for k in range (0,N):
                if Sx[k]==Sy[j] and Sx[k]==0:
                    psi=0   
                else:
                    psi=arccos(Sx[k]/(sqrt((Sx[k]**2)+(Sy[j]**2))))
                temp=0
                if a_or_r=="r":
                    temp=(1/(2*pi))*(cos(phi[i]-psi)+1)
                if a_or_r=="a":
                    temp=(1/(2*pi))*(-cos(phi[i]-psi)+1)
                Kernel_o_func[i,j,k]=temp
    return Kernel_o_func

