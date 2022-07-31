from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy.linalg import *
from compute_2D_FFT import *
from compute_2D_IFFT import *
from Matrix_A import *
from Matrix_B import *
from K_w import *
from S_u_combine import *
from Turning_Function import *
from newtimestep_M import *


N=64 #grid number 
sigma=pi/8
gamma=1 #speed
dt=0.01 #time step dt
tmax=3 #maximum time
q_r=0.5 #strength of repulsive force
q_a=2 #strength of attractive force
time_step=int(tmax/dt)
animation_use=array(zeros((time_step,N,N,N)))

'''
other matrices that does not change over time
'''
sign_matrix=sign_M(N)
K_w_a=K_w_g(N,"a")
K_w_r=K_w_g(N,"r")
K_w_fft_a=K_w_g_fft(N,K_w_a)
K_w_fft_r=K_w_g_fft(N,K_w_r)
A=matrixA(N,gamma,dt)
B=matrixB(N,gamma,dt)


print("Setup finished")
 
'''
Started loop over time
produce animation use matrix
that records u in different time
'''

for i in range (0, time_step):
    T_func=Turning(q_r, sign_matrix,K_w_fft_r,u, N)+Turning(q_a, sign_matrix,K_w_fft_a,u, N) #turning function
    Non_linear_func=Non_linear(N, T_func, u,3) #Non_linear term with/without saturation function
    u_new=newtimestep_matrix(u,Non_linear_func, gamma, dt,N, sign_matrix,A,B) #Calculate the new timestep u
    u=real(u_new) #update u
    animation_use[i]=u #record u
    print('range of u is [%e,%e], the range of T is [%e,%e].' % (amin(u),amax(u),amin(T_func),amax(T_func)))
print("time update finished")

'''
sum over then angle phi
'''
for i in range (0,time_step):
    for k in range (0,N):
        for l in range (0,N):
            sum_phi[i,k,l]=(2*pi/N)*sum(animation_use[i,:,k,l]) #normalized by 2*pi/N
save("u_file_name", animation_use)
save("sum_phi_file_name", sum_phi)
print("all data saved")
