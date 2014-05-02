from pyFi.quotes import api
from pyFi.greeks import volatility
from pyFi.methods import binomial as bi
from pyFi.methods import mc
from pyFi.methods import fd
from pyFi.chain import chain as ch

# This section gets our constants
exercise = 5.
print 'Exercise: $'+str(exercise)
spot = 4.
print 'Price:    $'+ str(spot)

sigma = 0.3
print 'Sigma:     '+ str(float(sigma)*100.0)+'%'

rfint = 0.04
print 'Rate:      '+str(float(rfint)*100.0)+'%'

T     = 1.
print 'Time:      '+str(T)

print '\n--------------------------'
#This section tests our solvers
solver = bi.binomial_euro(S=spot,E=exercise,r=rfint,M=60000,sigma=sigma,method='higham',T=T,opt='put')
print 'Binomial Method:    $'+str(solver.solve())

solver2 = mc.mcfast_euro(S=spot,E=exercise,r=rfint,M=100000000,sigma=sigma,T=T,opt='put')
print 'Monte Carlo Method: $'+str(solver2.solve())

solver3 = fd.BS_fd_explicit(S=spot,E=exercise,r=rfint,M=120,k=9000,L=exercise,sigma=sigma,T=T,opt='put')
print 'Finite Difference:  $'+str(float(solver3.solve()))
