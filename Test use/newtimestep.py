#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA

def newtimestep(u, Turning, Lambda, Moving, dt,N):
    ##u: latest value of u
    ##Turning,Lambda,Moving: these functions are 
    ##all calculated by using the latest valye of u
    sample_2D=zeros((N,N))
    sample_3D=[]
    for i in range (0,N):
        sample_3D.append(sample_2D)
    sample_3D=array(sample_3D)
    u_new=sample_3D
    du=sample_3D
    du=Turning-Lambda+Moving
    u_new=u+(du*dt)
    return u_new

#def newtimestepRK2(uq_hat, ur_hat, up_hat, F, gamma, dt, plus, minus):
    ###F=-G+G*2pi/n*w_hat
    ###plus:ik+l minus:ik-l
    ###uq_hat: current angle, ur_hat: next angle, up_hat_prev angle
    #k1=-gamma/2*(plus*ur_hat+minus*up_hat)+F*uq_hat
    #k2=-gamma/2*(plus*ur_hat+minus*up_hat)+F*(uq_hat+0.5*k1*dt)
    #uqnew_hat=uq_hat+k2*dt
    #return uqnew_hat
    
