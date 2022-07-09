from numpy import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (k and l)
## where [I-dtM]^(-1)

##Make sure k and l are integers but not values
def diagonal_matrixA(N,k,l,gamma,dt):
    w_number_use=int(N/2)
    w_number= arange(-w_number_use,w_number_use)+1
    w_number[-1]=0
    kk=w_number[k]
    ll=w_number[l]    
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    result=ones(N).astype(complex)
    for i in range (0,N):
        theta=phi[i]
        const=exp(-1j*gamma*(ll*cos(theta)+kk*sin(theta))*dt)
        #const=exp(-1j*gamma*(k*cos(theta)+l*(sin(theta)))*dt)
        result[i]=const
    return result
                
    
                
    
    
    
    