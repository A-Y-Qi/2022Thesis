from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from matplotlib import *



##data us a 4D array with the format (time, angle, x, y)
def Arrow_max(data,N):
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    timestep_number=len(data)
    sample_2D=zeros((N, N))
    sample_3D=[]
    for i in range (0,timestep_number):
        sample_3D.append(sample_2D)
    sample_3D=array(sample_3D)
    Ex=array(sample_3D)
    Ey=array(sample_3D)
    for i in range (0,timestep_number):
        for k in range (0,N):
            for l in range (0,N):
                temp_array=data[i,:,k,l]
                temp=max(temp_array)
                index_num=list(temp_array).index(temp)
                cos_x=cos(phi[index_num])
                sin_y=sin(phi[index_num])
                Ex[i,k,l]=cos_x
                Ey[i,k,l]=sin_y
    return Ex,Ey
    




def Arrow_avg(data,N):
    data=array(data)
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    timestep_number=len(data)
    sample_2D=zeros((N, N))
    sample_3D=[]
    for i in range (0,timestep_number):
        sample_3D.append(sample_2D)
    sample_3D=array(sample_3D)
    sample_3D_2=array(sample_3D)
    sample_3D_3=array(sample_3D)
    sample_3D_4=array(sample_3D)
    d_xy=sample_3D_2
    Ex=sample_3D_3
    Ey=sample_3D_4    
    for i in range (0,timestep_number):
        for k in range (0,N):
            for l in range (0,N):
                xy=data[i,:,k,l]
                xy_fft=fft.fft(xy)
                d_xy[i,k,l]=real(xy_fft[0])
    d_xy=array(d_xy)
    rho=array(data)
    comb=array(data)
    for i in range (0,timestep_number):
        for k in range (0,N):
            for l in range (0,N):
                density=d_xy[i,k,l]
                if density==0:
                    rho[i,:,k,l]=zeros(N)
                    print(i,k,l)
                else:
                    rho[i,:,k,l]=data[i,:,k,l]/density
    for j in range (0,N):
        angle=phi[j]
        comb[:,j,:,:]=angle*rho[:,j,:,:]
    for i in range (0,timestep_number):
        for k in range (0,N):
            for l in range (0,N):
                temp_array=comb[i,:,k,l]
                temp_fft=fft.fft(temp_array)
                temp_int=real(temp_fft[0])
                Ex[i,k,l]=cos(temp_int)
                Ey[i,k,l]=sin(temp_int)  
    return Ex,Ey
