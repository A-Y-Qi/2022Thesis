from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA



def fft_3D(matrixin3D, N):
    ##set up a new 3d matrix
    matrix_fft=matrixin3D
    ##start to find the dft of each rows first
    #take out each 2D matrix first
    for i in range (0, N):
        matrix_2D=matrixin3D[i]
        FFT_2D=matrix_2D
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for j in range(0, N):
            row=matrix_2D[j]
            FFT_row=fft.fftshift(fft.fft(row))
            ##Add normalization (2pi period)
            FFT_2D[j]=(2*pi/N)*FFT_row
        #take out and operate on the columns next
        for j in range(0, N):
            column=FFT_2D[:,j]
            FFT_column=fft.fftshift(fft.fft(column))
            ##Add normalization (2pi period)
            FFT_2D[:,j]=(2*pi/N)*FFT_column    
        matrix_fft[i]=FFT_2D
    #Now the matrix_fft has applied FFT twice
    #Take out the matrix_fft to apply for the third time
    for i in range (0,N):
        for j in range(0,N):
            height=matrix_fft[:,i,j]
            FFT_height=fft.fftshift(fft.fft(height))
            ##Add normalization (2pi period)
            matrix_fft[:,i,j]=(2*pi/N)*FFT_height
    return matrix_fft



            

