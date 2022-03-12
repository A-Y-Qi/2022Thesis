from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_4D_FFT import *

def Fourier_Turning(q,attract_kernel,turning_fun, u, N):
    a_kenel_fft=fft_3D(array(attract_kernel),N)
    turn_fft=fft_4D(array(turning_fun),N)
    fft_u=fft_3D(u,N)
    wave_0_u=fft_u[int(N/2)]
    for i in range (0,N):
        fft_u[i]=wave_0_u
    ##create 4D array
    temp=[]      
    for i in range (0,N):
        temp.append([])
        for j in range (0,N):
            temp[i].append([])
            for k in range(0,N):
                temp[i][j].append([])
                for l in range(0,N):
                    temp[i][j][k].append(1)
    Fourier_turning=array(temp).astype(complex)
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            #wavenumber for x
            for k in range (0, N):
                #wavenumber for y
                for l in range (0, N):
                    plus=mod(i+j,N)
                    temp2=q*a_kenel_fft[j,k,l]*turn_fft[i,plus,k,l]*fft_u[int(N/2),k,l]
                    Fourier_turning[i,j,k,l]=temp2
    return Fourier_turning
