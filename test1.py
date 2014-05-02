from pyFi.quotes import api
from pyFi.greeks import volatility
from pyFi.methods import binomial as bi
from pyFi.methods import mc
from pyFi.methods import fd
from pyFi.chain import chain as ch

# This section gets our constants
symbol = 'AAPL'
exercise = 590
print 'Exercise: $'+str(exercise)

spot = api.get_quote_google(symbol)
print symbol+ ' Price: $'+ str(spot)

sigma = volatility.hist_vol(symbol,30)
print symbol+' Sigma: '+ str(float(sigma)*100.0)+'%'

rfint = api.get_risk_free()
print 'Risk Free rate: '+str(float(rfint)*100.0)+'%'
print '\n--------------------------'

#This section tests our solvers
solver = bi.binomial_euro(S=spot,E=exercise,r=rfint,M=40000,sigma=sigma,method='higham',T=2.0/365.0,opt='put')
print 'Binomial Method:    $'+str(solver.solve())

solver2 = mc.mcfast_euro(S=spot,E=exercise,r=rfint,M=4000000,sigma=sigma,T=2.0/365.0,opt='put')
print 'Monte Carlo Method: $'+str(solver2.solve())

solver3 = fd.BS_fd_explicit(S=spot,E=exercise,r=rfint,M=1201,k=1700,L=2*exercise,sigma=sigma,T=2.0/365.0,opt='put')
print 'Finite Difference:  $'+str(float(solver3.solve()))

#This section produces an option chain
chain = ch.option_chain(spot,rfint,sigma,T=2.0/365,range=(520,5,640),method='bi')
chain.generate()
chain.export_csv('test.csv')
print '\n--------------------------\nOpen test.csv for option chain.\n'