from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA



def ifft_3D(FFT_matrixin3D, N):
    ##set up a new 3d matrix
    matrix_ifft=FFT_matrixin3D
    ##start to find the dft of each rows first
    #take out each 2D matrix first
    for i in range (0, N):
        FFT_2D=FFT_matrixin3D[i]
        matrix_2D=FFT_2D
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for j in range(0, N):
            row=FFT_2D[j]
            matrix_row=fft.ifft(fft.ifftshift(row))
            ##Add normalization (2pi period)
            matrix_2D[j]=(2*pi/N)*matrix_row
        #take out and operate on the columns next
        for j in range(0, N):
            column=matrix_2D[:,j]
            matrix_column=fft.ifft(fft.ifftshift(column))
            ##Add normalization (2pi period)
            matrix_2D[:,j]=(2*pi/N)*matrix_column    
        matrix_ifft[i]=matrix_2D
    #Now the matrix_fft has applied FFT twice
    #Take out the matrix_fft to apply for the third time
    for i in range (0,N):
        for j in range(0,N):
            height=matrix_ifft[:,i,j]
            matrix_height=fft.ifft(fft.ifftshift(height))
            ##Add normalization (2pi period)
            matrix_ifft[:,i,j]=(2*pi/N)*matrix_height
    return matrix_ifft



            

