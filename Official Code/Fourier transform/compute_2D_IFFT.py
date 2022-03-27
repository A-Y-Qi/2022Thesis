from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA



def ifft_2D(matrixin2D, N):
        ##set up a new 3d matrix
        IFFT_temp=matrixin2D
        IFFT_2D=zeros((N, N))
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for j in range(0, N):
                row=matrixin2D[j]
                IFFT_row=zeros(N)
                IFFT_row=fft.ifft(fft.ifftshift(row))
                #IFFT_row=fft.ifft(row)
                ##Add normalization (2pi period)
                IFFT_row_real=real(IFFT_row)
                IFFT_2D[j]=(N/(2*pi))*IFFT_row_real
                #take out and operate on the columns next
        for j in range(0, N):
                column=IFFT_2D[:,j]
                IFFT_column=zeros(N)
                IFFT_column=fft.ifft(fft.ifftshift(column))
                #IFFT_column=fft.ifft(column)
                ##Add normalization (2pi period)
                IFFT_column_real=real(IFFT_column)
                IFFT_2D[:,j]=(N/(2*pi))*IFFT_column_real
                #Now the matrix_fft has applied FFT twice
        return IFFT_2D



            

