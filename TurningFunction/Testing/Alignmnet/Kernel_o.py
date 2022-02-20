from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special

def Kernel_o(N):
    #define theta and phi
    theta = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #Generate a 2*2 matrix
    Kernel_o_func = zeros((N,N))
    for i in range (0,N):
        for j in range(0,N):
            Kernel_o_func[i][j]=(1/(2*pi))*(cos(phi[i]-theta[j])+1)
    return Kernel_o_func