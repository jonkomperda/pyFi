from numpy import *

E = 4
sigma = 0.3
r = 0.03
T = 1.
L = 10.
M = 11
k = 10
dt = T/k
h = L/(M-1)

n = array([x for x in range(0,M)])

a = 0.5*dt*(sigma**2 * n**2 - r*n)
b = 1.0 - dt*(sigma**2 * n**2 + r)
c = 0.5*dt*(sigma**2*n**2 + r*n)

Mat = diag(b) + diag(c[0:M-1],1) + diag(a[1:M],-1)
Mat[0,:] = 0
Mat[0,0] = 1
Mat[M-1,:] = 0
Mat[M-1,M-1] = 1 
print Mat

S = h*n
print S
Uold = reshape(maximum(E-S,0),(M,1))
print Uold

for i in range(k):
    U = dot(Mat,Uold)
    U[0,0] = E*exp(-r*dt*(i+1))
    print U
    Uold = U




