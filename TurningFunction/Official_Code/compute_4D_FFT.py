from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from compute_3D_FFT import *


def fft_4D(matrixin4D, N):
    ##set up a new 4d matrix
    matrix4D_fft=array(array(matrixin4D).astype(complex))
    result=matrix4D_fft
    ##start to find the dft of each rows first
    #take out each 3D matrix first
    for l in range (0, N):
        matrix_3D=array(matrixin4D[l])
        temp=array(fft_3D(matrix_3D,N))
        matrix4D_fft[l]=temp
    for i in range (0,N):
        for j in range (0,N):
            for k in range(0,N):
                height_4D=matrix4D_fft[:,i,j,k]
                FFT_height_4D=array(fft.fftshift(fft.fft(height_4D)))
                result[:,i,j,k]=((2*pi)/N)*FFT_height_4D
    return result
