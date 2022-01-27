from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA


def fft_3D(matrixin3D, N):
    ##set up a new 3d matrix
    matrix_fft=matrixin3D
    ##start to find the dft of each rows first
    #take out each 2D matrix first
    for i in range (0, N-1):
        matrix_2D=matrixin3D[i]
        FFT_2D=matrix_2D
        ##then apply FFT to each row and column
        #take out and operate on the rows first
        for j in range(0, N-1):
            row=matrix_2D[j]
            FFT_row=fft.fft(row)
            ##Add normalization
            FFT_2D[j]=1/(sqrt(N))*FFT_row
        #take out and operate on the columns next
        for j in range(0, N-1):
            column=FFT_2D[:,j]
            FFT_column=fft.fft(column)
            ##Add normalization
            FFT_2D[:,j]=1/(sqrt(N))*FFT_column    
        matrix_fft[i]=FFT_2D
    #Now the matrix_fft has applied FFT twice
    #Take out the matrix_fft to apply for the third time
    for i in range (0,N-1):
        for j in range(0,N-1):
            height=matrix_fft[:,i,j]
            FFT_height=fft.fft(height)
            ##Add normalization
            matrix_fft[:,i,j]=1/(sqrt(N))*FFT_height
    return matrix_fft

##The following part is for the testing part
##Use the FFT paper code for testing
###Normalization for the fft.fft is 1/(sqrt(N))

##Construct a testing matrix
##2cos(x)*2sin(y)*(2cos(z))^2

#set up a 3D matrix of zero
N=10
sample_2D=zeros((N, N),dtype=complex)
sample_3D=[]
for i in range (0,N-1):
    sample_3D.append(sample_2D)
sample_3D=array(sample_3D)

##Construct a testing matrix
##2cos(x)*2sin(y)*(2cos(z))^2
h=2*pi/N
for i in range (0,N-1):
    for j in range (0,N-1):
        for k in range(0,N-1):
            sample_3D[i][j][k]=2*cos(i*h)*2*sin(j*h)*(2*cos(k*h))**2

            

