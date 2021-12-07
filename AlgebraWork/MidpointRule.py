from numpy import *

def midpoint(f, a, b, n):
    h = float(b-a)/n
    result = 0
    for i in range(n):
        result += f((a + h/2.0) + i*h)
    result *= h
    return result

g = lambda y: exp(-y**2)
a=0
b=1
n=10

k=midpoint(g,0,1,10)