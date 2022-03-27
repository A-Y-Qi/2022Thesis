from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA



def fft_2D(matrixin2D, N):
        ##set up a new 3d matrix
        FFT_temp=matrixin2D
        FFT_2D=zeros((N, N),dtype=complex)
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for j in range(0, N):
                row=matrixin2D[j]
                FFT_row=zeros(N, dtype=complex)
                FFT_row=fft.fftshift(fft.fft(row))
                ##Add normalization (2pi period)
                FFT_2D[j]=(2*pi/N)*FFT_row
                #take out and operate on the columns next
        for j in range(0, N):
                column=FFT_2D[:,j]
                FFT_row=zeros(N, dtype=complex)
                FFT_column=fft.fftshift(fft.fft(column))
                ##Add normalization (2pi period)
                FFT_2D[:,j]=(2*pi/N)*FFT_column
                #Now the matrix_fft has applied FFT twice
        return FFT_2D



            

