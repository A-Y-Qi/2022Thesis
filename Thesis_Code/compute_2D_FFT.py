from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA

def sign_M(N):
    sign_matrix=zeros((N, N))
    for i in range (N):
        for j in range (N):
            sign_matrix[i,j]=(-1)**(i+j)
    return sign_matrix


#No FFTshift Applied
def fft_2D_3DM(matrixin3D, sign_matrix,N):
        ##set up a new 3d matrix
        matrix_fft= empty_like(matrixin3D,dtype=cdouble)
        for i in range (0,N):
                matrix_fft[i,:,:]=(((1/float(N)))**2)*fft.fft2(matrixin3D[i,:,:])*sign_matrix
        return matrix_fft

#No FFTshift Applied
def fft_2D_4DM(matrixin4D,sign_matrix, N):
        ##set up a new 3d matrix
        matrix_fft=array(matrixin4D).astype(complex)
        for i in range (0,N):
                for k in range (0,N):
                        #matrix_fft[i,k,:,:]=(((2*pi/N))**2)*fft.fftshift(fft.fft2(matrixin4D[i,k,:,:]))
                        matrix_fft[i,k,:,:]=(((1/N))**2)*fft.fft2(matrixin4D[i,k,:,:])*sign_matrix
                
        return matrix_fft