from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA


def sign_M(N):
    sample_2D=zeros((N, N))
    sign_matrix=copy(sample_2D)
    w_number_sign=fft.ifftshift((-1)**(arange(-N/2.0,N/2.0)+1))
    for i in range (0,N):
        sign_matrix[:,i]=w_number_sign[i]*w_number_sign
    return sign_matrix

def fft_2D(matrixin2D, N):
        ##set up a new 3d matrix
        FFT_temp=matrixin2D
        FFT_2D=zeros((N, N),dtype=complex)
        zero_wnum_index=int((N/2)-1)
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for j in range(0, N):
                row=matrixin2D[j]
                FFT_row=zeros(N, dtype=complex)
                FFT_row=fft.fftshift(fft.fft(row))
                #FFT_row[-1]=FFT_row[zero_wnum_index]
                ##Add normalization (2pi period)
                FFT_2D[j]=(2*pi/N)*FFT_row
                #take out and operate on the columns next
        for j in range(0, N):
                column=FFT_2D[:,j]
                FFT_row=zeros(N, dtype=complex)
                FFT_column=fft.fftshift(fft.fft(column))
                #FFT_column[-1]=FFT_column[zero_wnum_index]
                ##Add normalization (2pi period)
                FFT_2D[:,j]=(2*pi/N)*FFT_column
                #Now the matrix_fft has applied FFT twice
        return FFT_2D

def fft_2D_3DM(matrixin3D, sign_matrix, N):
        ##set up a new 3d matrix
        matrix_fft=array(matrixin3D).astype(complex)
        for i in range (0,N):
                matrix_fft[i,:,:]=(((2*pi/N))**2)*fft.fftshift(fft.fft2(matrixin3D[i,:,:]))*sign_matrix
        return matrix_fft

#No Normalization Applied
def fft_2D_4DM(matrixin4D,sign_matrix, N):
        ##set up a new 3d matrix
        matrix_fft=array(matrixin4D).astype(complex)
        for i in range (0,N):
                for k in range (0,N):
                        #matrix_fft[i,k,:,:]=(((2*pi/N))**2)*fft.fftshift(fft.fft2(matrixin4D[i,k,:,:]))
                        matrix_fft[i,k,:,:]=(((2*pi/N))**2)*fft.fft2(matrixin4D[i,k,:,:])*sign_matrix
                
        return matrix_fft
