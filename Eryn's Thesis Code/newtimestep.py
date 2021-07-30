#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA

def newtimestep(uq_hat, ur_hat, up_hat, F, gamma, dt, plus, minus):
    ##F=-G+G*2pi/n*w_hat
    ##plus:ik+l minus:ik-l
    ##uq_hat: current angle, ur_hat: next angle, up_hat_prev angle
    uqnew_hat=(-gamma/2*plus*ur_hat-gamma/2*minus*up_hat+F*uq_hat)*dt+uq_hat
    return uqnew_hat

def newtimestepRK2(uq_hat, ur_hat, up_hat, F, gamma, dt, plus, minus):
    ##F=-G+G*2pi/n*w_hat
    ##plus:ik+l minus:ik-l
    ##uq_hat: current angle, ur_hat: next angle, up_hat_prev angle
    k1=-gamma/2*(plus*ur_hat+minus*up_hat)+F*uq_hat
    k2=-gamma/2*(plus*ur_hat+minus*up_hat)+F*(uq_hat+0.5*k1*dt)
    uqnew_hat=uq_hat+k2*dt
    return uqnew_hat
    
