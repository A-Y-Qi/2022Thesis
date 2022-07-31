from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_2D_FFT import *
from scipy import signal as sig

#q: a constant
#K_w_fft: the K_w function with its fft form over X, ie K_w_fft=fft_2D_4DM(K_w_func,N)
#turning_fun: 2D for alignment and 4D for attraction/repulsion
#u: a 3D function

def Turning(q, sign_matrix, K_w_fft,u, N):
    ##create 4D array
    temp_2D=zeros((N,N))
    temp_3D=zeros((N,N,N))
    Fourier_turning=array(zeros((N,N,N,N))).astype(complex)
    turning=array(zeros((N,N,N,N)))
    ##Integral over theta
    sum_phi=copy(temp_2D).astype(complex)
    for k in range (0,N):
        for l in range (0,N):
            sum_phi[k,l]=(2*pi/N)*real(fft.fft(u[:,k,l])[0])
    sum_phi_fft=fft.fft2(sum_phi)*sign_matrix*((1/N)**2)
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            Fourier_turning[i,j,:,:]=K_w_fft[i,j,:,:]*sum_phi_fft*sign_matrix
            turning[i,j,:,:]=q*((N/(1))**2)*real(fft.ifft2(Fourier_turning[i,j,:,:]))
    return turning