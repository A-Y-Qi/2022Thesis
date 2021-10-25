from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from grid_setup import *
from fouriergrid_setup import *
from newtimestep import *
from angle_list import *
from plot_animation import *

b=5
test_matrix = random.randint(10, size=(b,b))
print("test_matrix",test_matrix)
matrix_fft=fft.fft2(test_matrix)
row_fft = zeros((b,b), dtype=complex)
column_fft=zeros((b,b), dtype=complex)

#sample original matrix
sample_3D=[]
for i in range (0,b):
    w=random.randint(10, size=(b,b))
    sample_3D.append(w)
sample_3D=array(sample_3D)
print("sample_3D=",sample_3D)

#new temp matrix
z_fft=[]      
for i in range (0,b):
    z_fft.append([])
    for j in range (0,b):
        z_fft[i].append([])
        for k in range(0,b):
            z_fft[i][j].append(0)
print("z_fft=",z_fft)
    
    
for k in range(0,b):
    test_matrix=sample_3D[k]
    for i in range(0,b):
        row_fft_i=fft.fft(test_matrix[i])
        row_fft[i]=row_fft_i
    for i in range(0,b):
        column_fft_i=fft.fft(row_fft[:,i])
        column_fft[:,i]=column_fft_i    
    z_fft[k]=column_fft
z_fft=array(z_fft)

for i in range(0,b):
    for j in range(0,b):
        z_fft_k=fft.fft(z_fft[:,i,j])
        z_fft[:,i,j]=array(z_fft_k)
    

sample_3D[:,0,0]