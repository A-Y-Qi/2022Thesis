from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA



def fft_3D(matrixin3D, N):
    ##set up a new 3d matrix
    matrix_fft=array(matrixin3D).astype(complex)
    ##start to find the dft of each rows first
    #take out each 2D matrix first
    for i in range (0, N):
        matrix_2D=array(matrixin3D[i])
        FFT_2D=matrix_2D.astype(complex)
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


####test on FFT
N=32
m = (2*pi/N)*arange(-N/2.0,N/2.0)
sample_2D=zeros((N, N),dtype=complex)
sample_3D=[]
for i in range (0,N):
    sample_3D.append(sample_2D)
sample_3D=array(sample_3D)
sample_3D_2=sample_3D
sample_3D_3=sample_3D

cos_5=cos(5*m)
cos_4=cos(4*m)
cos_2=cos(2*m)

##Construct a testing matrix
##2cos(5x)*2sin(4y)*(2cos(z))^2
h=2*pi/N
for i in range (0,N):
    for j in range (0,N):
        for k in range(0,N):
            temp=cos_5[i]*cos_4[j]*cos_2[k]
            sample_3D[i][j][k]=temp
wc=fft_3D(sample_3D, N)

            

