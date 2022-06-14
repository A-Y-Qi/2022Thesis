from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_2D_FFT import *
from compute_2D_IFFT import *
from Kw_combined import *
from scipy import signal as sig

#q: a constant
#K_w_fft: the K_w function with its fft form over X, ie K_w_fft=fft_2D_4DM(K_w_func,N)
#turning_fun: 2D for alignment and 4D for attraction/repulsion
#u: a 3D function

def Turning(q, sign_matrix, K_w_fft,u, N):
    ##create 4D array
    temp_2D=zeros((N,N))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    Fourier_turning=[]
    turning=[]
    for i in range (0,N):
        array(Fourier_turning.append(temp_3D))
        array(turning.append(temp_3D))
    Fourier_turning=array(Fourier_turning).astype(complex)
    turning=array(turning).astype("float64")
    ##Integral over theta
    sum_phi=copy(temp_2D).astype(complex)
    for k in range (0,N):
        for l in range (0,N):
            sum_phi[k,l]=(2*pi/N)*real(fft.fft(u[:,k,l])[0])
    sum_phi_fft=fft.fft2(sum_phi)*sign_matrix
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            Fourier_turning[i,j,:,:]=K_w_fft[i,j,:,:]*sign_matrix*sum_phi_fft
            turning[i,j,:,:]=real(fft.ifft2(Fourier_turning[i,j,:,:]))
    return turning
            
    
                    
            
         
                   
                    
                    
