from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy.linalg import *
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_3D_IFFT import *
from compute_2D_IFFT import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (k and l)
## where [I-dtM]^(-1)

##Make sure k and l are integers but not values
def diagonal_matrixA(N,k,l,gamma,dt):
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    result=ones(N).astype(complex)
    for i in range (0,N):
        theta=phi[i]
        const=exp(-1j*gamma*(k*cos(theta)+l*sin(theta))*dt)
        #const=exp(-1j*gamma*(k*cos(theta)+l*(sin(theta)))*dt)
        result[i]=const
    return result
