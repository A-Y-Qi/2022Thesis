from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy.linalg import *
import random
#from compute_3D_FFT import *
#from compute_2D_FFT import *
#from compute_3D_IFFT import *
#from compute_2D_IFFT import *
#from Matrix_A import *
#from Matrix_B import *
#from Kw_combined_debug import *
#from S_u_combine_debug import *
#from Debug_use_Turning_Function import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (k and l)
## where [I-dtM]^(-1)

"""
This document is a test document which record the test between
my code and lennaert's code
meanwhile I may use this document to generate initial condition 
for comparing
"""

'''
Start with our first test: single Gaussian bump
'''
#x0 = 0.4
#y0 = -0.2
#sigx = 0.7
#sigy = 0.9
#phi0 = pi*0.9
#sigphi = 0.5
#mass=10000
#for i in range(0,N):
    #for k in range(0,N):
        #for l in range(0,N):
            #u[i,k,l] = mass*g(position[l]-x0,sigx,N)*g(position[k]-y0,sigy,N)*g(position[i]-phi0,sigphi,N)

#save("1_bump_09pi_Lennaert", u)
u=load("1_bump_09pi_Lennaert.npy")
'''
compare the K_w function first:
***My original weight function is delta
Now I changed it to periodic gaussian
which is the same as lennaert's to compare
'''
my_kw_fft_a=load("K_w_fft_64_attraction_Lennaert.npy")
my_kw_fft_r=load("K_w_fft_64_attraction_Lennaert.npy")
len_kw_fft_a=load("Qa_fft_Lennaert.npy")
len_kw_fft_r=load("Qr_fft_Lennaert.npy")

###############################################
#Difference:
#linalg.norm(my_kw_fft_a)
#9.726552472178765
#linalg.norm(len_kw_fft_a)
#0.2464577778009756
#linalg.norm(my_kw_fft_r)
#9.726552472178765
#linalg.norm(len_kw_fft_r)
#2.019765577493544
###############################################


'''
parameters:
1. Attraction:
qa=1, ma=0.5, da=0.5, ka=1
mass=10000
N=64 X 64 X 64
st=0.05
'''
