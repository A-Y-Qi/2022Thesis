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
def linear_matrix(N,k,l,gamma,dt):
    w_number=arange(-N/2.0,N/2.0)
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    identity_matrix=identity(N)
    kk=w_number[k]
    ll=w_number[l]
    result=zeros((N,N),dtype=complex)
    M=zeros((N,N),dtype=complex)
    temp=zeros((N,N),dtype=complex)
    column=zeros(N)*1j
    row=zeros(N)*1j
    column[1]=(-gamma/2)*(1j*kk+ll)
    column[-1]=(-gamma/2)*(1j*kk-ll)
    row[1]=(-gamma/2)*(1j*kk-ll)
    row[-1]=(-gamma/2)*(1j*kk+ll)    
    M=toeplitz(column, row)
    M=dt*M
    temp=identity_matrix-M
    result=inv(temp)
    return result
                
    
                
    
    
    
    