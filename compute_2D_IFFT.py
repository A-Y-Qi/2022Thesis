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

def ifft_2D(matrixin2D, N):
        ##set up a new 3d matrix
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
        
        
        
def ifft_2D_3DM(matrixin3D, sign_matrix, N):
        ##set up a new 3d matrix
        IFFT_3D=array(real(matrixin3D))
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for i in range (0,N):
                IFFT_3D[i,:,:]=(((N/(2*pi)))**2)*real(fft.ifft2(fft.ifftshift(matrixin3D[i,:,:]*sign_matrix)))
        return IFFT_3D

##No fftshift applied
def ifft_2D_4DM(matrixin4D,sign_matrix, N):
        ##set up a new 3d matrix
        IFFT_4D=array(real(matrixin4D))
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for i in range (0,N):
                for k in range (0,N):
                        IFFT_4D[i,k,:,:]=(((N/(2*pi)))**2)*real(fft.ifft2(matrixin4D[i,k,:,:])*sign_matrix)
        return IFFT_4D



            

