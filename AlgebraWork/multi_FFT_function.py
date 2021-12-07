#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA


def fft_3D(matrixin3D,size):
    row_fft = zeros((size,size), dtype=complex)
    column_fft=zeros((size,size), dtype=complex)
    z_fft=[]      
    for i in range (0,size):
        z_fft.append([])
        for j in range (0,size):
            z_fft[i].append([])
            for k in range(0,size):
                z_fft[i][j].append(1) 
    #fft loop
    for k in range(0,size):
        test_matrix=matrixin3D[k]
        for i in range(0,size):
            row_fft_i=fft.fft(test_matrix[i])
            row_fft[i]=row_fft_i
        for i in range(0,size):
            column_fft_i=fft.fft(row_fft[:,i])
            column_fft[:,i]=column_fft_i    
        z_fft[k]=column_fft
    z_fft=array(z_fft)
    for i in range(0,size):
        for j in range(0,size):
            z_fft_k=fft.fft(z_fft[:,i,j])
            z_fft[:,i,j]=array(z_fft_k)
    return z_fft

def fft_2D(matrixin2D,size):
    row_fft = zeros((size,size), dtype=complex)
    column_fft=zeros((size,size), dtype=complex)    
    for i in range(0,size):
        row_fft_i=fft.fft(matrixin2D[i])
        row_fft[i]=row_fft_i
    for i in range(0,size):
        column_fft_i=fft.fft(row_fft[:,i])
        column_fft[:,i]=column_fft_i 
    return column_fft