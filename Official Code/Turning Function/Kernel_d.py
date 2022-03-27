from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special

#range over 2pi
def Kernel_d(N,m,d):
    #define x and y
    x = (2*pi/N)*arange(-N/2.0,N/2.0)
    y = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define constant A
    A=m*(m*exp(-(d**2)/(m**2))+d+d*special.erf(d/m))
    #Generate a 2*2 matrix
    Kernel_d_func = zeros((N,N))
    for i in range (0,N):
        for j in range(0,N):
            Kernel_d_func[i,j]=(2/A)*exp((-(sqrt(x[i]**2+y[j]**2)-d)**2)/(m**2))
    return Kernel_d_func

    
