from numpy import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (k and l)
## where [I-dtM]^(-1)

##Make sure k and l are integers for wavenumbers but not values
def diagonal_matrixA(N,k,l,gamma,dt):
    w_number= concatenate((linspace(0,N/2,int(N/2+1)),linspace(-N/2+1,-1,int(N/2-1))))
    w_number[int(N/2)] = 0
    kk=w_number[k]
    ll=w_number[l]    
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    result=ones(N).astype(complex)
    for i in range (0,N):
        theta=phi[i]
        const=exp(-1j*gamma*(ll*cos(theta)+kk*sin(theta))*dt)
        result[i]=const
    return result

def matrixA(N,gamma,dt):
    result=array(zeros((N,N,N))).astype(complex)
    w_number = concatenate((linspace(0,N/2,int(N/2+1)),linspace(-N/2+1,-1,int(N/2-1))))
    w_number[int(N/2)] = 0
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    for k in range (0,N):
        for l in range (0,N):
            M = - (w_number[l] * cos(phi) + w_number[k] * sin(phi))
            result[:,k,l]=exp(dt*gamma* 1J * M)
    return result
                
    
                
    
    
    
    