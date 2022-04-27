from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy.linalg import *
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_3D_IFFT import *
from compute_2D_IFFT import *
from Matrix_A import *
from Matrix_B import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (k and l)
## where [I-dtM]^(-1)

##Make sure k and l are integers but not values
def newtimestep_matrix(u, Non_lin, gamma,dt,N):
    #first to transfer u to u_hat with only i and j
    u_fft=fft_2D_3DM(u, N)
    Non_lin_fft=fft_2D_3DM(Non_lin, N)
    u_new_fft=copy(u_fft)
    u_new=copy(u)
    for k in range (0,N):
        for l in range (0,N):
            u_hat_phi=u_fft[:,k,l]
            Non_lin_hat_phi=Non_lin_fft[:,k,l]
            matrixA=diagonal_matrixA(N,k,l,gamma,dt)
            matrixB=diagonal_matrixB(N,k,l,gamma,dt)
            u_new_phi=copy(u_hat_phi)
            #u_new_phi=matmul(matrixA, u_hat_phi)+matmul(matrixB, Non_lin_hat_phi)
            u_new_phi=matrixA*u_hat_phi+matrixB*Non_lin_hat_phi
            u_new_fft[:,k,l]=u_new_phi
    u_new=real(ifft_2D_3DM(u_new_fft, N))
    return u_new
    
                
    
                
    
    
    
    