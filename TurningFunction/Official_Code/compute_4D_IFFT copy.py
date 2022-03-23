from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from compute_3D_FFT import *


def fft_4D(matrixin4D, N):
    ##set up a new 4d matrix
    matrix_fft=array(matrixin4D)
    result=matrix_fft
    ##start to find the dft of each rows first
    #take out each 3D matrix first
    for i in range (0, N):
        matrix_3D=array(matrixin4D[i])
        matrix_fft[i]=array(fft_3D(matrix_3D,N))
    for i in range (0,N):
        for j in range (0,N):
            for k in range(0,N):
                height=matrix_fft[:,i,j,k]
                IFFT_height=array(fft.fftshift(fft.fft(height)))
                result[:,i,j,k]=((2*pi)/N)*IFFT_height
    return result



            

