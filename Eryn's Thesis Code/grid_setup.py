from numpy import *


def grid_setup(N):
    if (N == 0):
        print ("unable to setup grid")
    else:
        ##setup spacial position grid
        x = (2/N)*arange(-N/2.0,N/2.0)
        y = (2/N)*arange(-N/2.0,N/2.0)        
        xx,yy=meshgrid(x,y)
        xx=xx.flatten()
        yy=yy.flatten()
        return xx, yy