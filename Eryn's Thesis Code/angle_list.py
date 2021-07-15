#input the libraries
from numpy import *

def angle_list(phi_num,u_list,N):
    #the number of the angle position 
    #The big 2D list
    #N: gird number for x and y axis
    k=len(u_list)-1
    density_list=[]
    for i in range(0, k+1):
        density_phi=u_list[i][phi_num]
        density_phi=density_phi.reshape((N, N))
        density_list.append(list(density_phi))
    return density_list
    