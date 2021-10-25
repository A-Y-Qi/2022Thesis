from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from grid_setup import *
from fouriergrid_setup import *
from newtimestep import *
from angle_list import *
from plot_animation import *

test_matrix = random.randint(10, size=(20,20))
matrix_fft=fft.fft2(test_matrix)
row_fft = zeros((20,20), dtype=complex)
column_fft=zeros((20,20), dtype=complex)

##
for i in range(0,20):
    row_fft_i=fft.fft(test_matrix[i])
    row_fft[i]=row_fft_i
    
for i in range(0,20):
    column_fft_i=fft.fft(row_fft[:,i])
    column_fft[:,i]=column_fft_i
    
K=matrix_fft-column_fft
print(K)
    
