from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_3D_FFT import *
from compute_2D_FFT import *

#q: a constant
#dist_kernel: distance kernel, a 2D array
#orent_kernel: orientation kernel, a 1D array for alignment, 3D array for repulsion and attraction
#turning_fun: 2D for alignment and 4D for attraction/repulsion
#u: a 3D function
def Fourier_Turning(q,dist_kernel, orent_kernel,turning_fun, u, N):
    d_kenel_fft=fft_2D(dist_kernel, N)
    o_kenel_fft=fft.fftshift(fft.fft(orent_kernel))
    turn_fft=fft_2D(turning_fun,N)
    fft_u=fft_3D(u,N)
    ##create 4D array
    temp=[]      
    for i in range (0,N):
        temp.append([])
        for j in range (0,N):
            temp[i].append([])
            for k in range(0,N):
                temp[i][j].append([])
                for l in range(0,N):
                    temp[i][j][k].append(0)
    Fourier_turning=temp
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            #wavenumber for x
            for k in range (0, N):
                #wavenumber for y
                for l in range (0, N):
                    plus=mod(k+l,N)
                    Fourier_turning[i][j][k][l]=q*d_kenel_fft[k][l]*o_kenel_fft[i]*turn_fft[i][j]*fft_u[plus][k][l]
    return Fourier_turning
