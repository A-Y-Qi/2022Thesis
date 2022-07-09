from numpy import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (l and k)
## where [I-dtM]^(-1)

##Make sure k and l are integers but not values
def diagonal_matrixB(N,k,l,gamma,dt):
    w_number_use=int(N/2)
    #w_number= fft.fftshift(arange(-w_number_use,w_number_use)+1)
    w_number= arange(-w_number_use,w_number_use)+1
    w_number[-1]=0
    kk=w_number[k]
    ll=w_number[l]     
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    identity_matrix=identity(N)
    result=zeros(N).astype(complex)
    for i in range (0,N):
        theta=phi[i]
        M=-1j*gamma*(ll*(cos(theta))+kk*(sin(theta)))
        const=1j*dt*exp(dt*M*0.5)*sinc(dt*0.5*M/pi)##note sinc in python has sin(pi*x)/x
        result[i]=const
    return result
                
    
                
    
    
    
    