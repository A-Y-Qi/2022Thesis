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
        
##No fftshift applied        
def ifft_2D_3DM(matrixin3D, sign_matrix, N):
        ##set up a new 3d matrix
        IFFT_3D=array(real(matrixin3D))
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for i in range (0,N):
                IFFT_3D[i,:,:]=(((N/(1)))**2)*real(fft.ifft2(matrixin3D[i,:,:]*sign_matrix))
        return IFFT_3D

##No fftshift applied
def ifft_2D_4DM(matrixin4D,sign_matrix, N):
        ##set up a new 3d matrix
        IFFT_4D=array(real(matrixin4D))
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for i in range (0,N):
                for k in range (0,N):
                        IFFT_4D[i,k,:,:]=(((N/(1)))**2)*real(fft.ifft2(matrixin4D[i,k,:,:]*sign_matrix))
        return IFFT_4D



            

