from numpy import *

def fouriergrid_setup(N):
    if (N == 0):
        print ("unable to setup grid")
    else:
        ##setup fourier spacial position grid
        k=hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
        l=hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
        kk,ll=meshgrid(k,l)
        kk=kk.flatten()
        ll=ll.flatten()
        plus=1j*kk+ll
        minus=1j*kk-ll
        return kk, ll, plus, minus