from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy.linalg import *
from compute_2D_FFT import *
from compute_2D_IFFT import *

'''
with known matrix A B 
'''
def newtimestep_matrix(u,Non_lin, gamma,dt,N,sign_matrix,A,B):
    #first to transfer u to u_hat with only i and j
    u_fft=fft_2D_3DM(u, sign_matrix,N)
    Non_lin_fft=fft_2D_3DM(Non_lin, sign_matrix, N)
    u_new_fft=copy(u_fft)
    u_new=copy(u)
    u_new_fft=A*u_fft+B*Non_lin_fft
    u_new=real(ifft_2D_3DM(u_new_fft, sign_matrix,N))
    return u_new


def newtimestep_matrix_linear(u, gamma,dt,N,sign_matrix):
    #first to transfer u to u_hat with only i and j
    u_fft=fft_2D_3DM(u, sign_matrix,N)
    u_new_fft=copy(u_fft)
    u_new=copy(u)
    for k in range (0,N):
        for l in range (0,N):
            u_hat_phi=u_fft[:,k,l]
            matrixA=diagonal_matrixA(N,k,l,gamma,dt)
            u_new_fft[:,k,l]=matrixA*u_hat_phi
    u_new=real(ifft_2D_3DM(u_new_fft, sign_matrix,N))
    return u_new

