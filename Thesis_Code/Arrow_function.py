from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from matplotlib import *
            
def Arrow_max1000(data,sum_phi,N):
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    timestep_number=len(data)
    sample_2D=zeros((timestep_number, 1000))
    Ex=array(sample_2D)
    Ey=array(sample_2D)
    x=array(sample_2D)
    y=array(sample_2D)    
    for i in range (0,timestep_number):
        ind = argpartition(sum_phi[i], -1000,axis=None)[-1000:]
        for index in range(1000):
            coord=ind[index]
            y[i,index]=phi[coord//N]
            x[i,index]=phi[coord%(N)]
            Ex[i,index]=cos(phi[argmax(data[i,:,coord//N,coord%(N)])])
            Ey[i,index]=sin(phi[argmax(data[i,:,coord//N,coord%(N)])])
    return Ex,Ey, x,y  

def Arrow_max2000(data,N):
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    timestep_number=len(data)
    sample_2D=zeros((timestep_number, 2000))
    Ex=array(sample_2D)
    Ey=array(sample_2D)
    x=array(sample_2D)
    y=array(sample_2D)    
    for i in range (0,timestep_number):
        ind = argpartition(data[i], -2000,axis=None)[-2000:]
        for index in range(2000):
            coord=ind[index]
            Ex[i,index]=cos(phi[coord//(N**2)])
            Ey[i,index]=sin(phi[coord//(N**2)])
            y[i,index]=phi[(coord%(N**2))//N]
            x[i,index]=phi[coord%(N)]
    return Ex,Ey, x,y  
