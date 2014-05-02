import sys
sys.path.append('../')
try:
    from ..methods import fd
    from ..methods import binomial as bi
    from ..methods import mc
except:
    import methods.fd as fd
    import methods.mc as mc
    import methods.binomial as bi

## Main class for the creation of the option chain. Must be called for initialization
class option_chain():
    ## Called upon initialization of the option chain
    # @param S = Spot Price of asset
    # @param r = Risk free interest rate
    # @param sigma = Volatility of the asset (\f$\sigma\f$)
    # @param T = Time to expiry
    # @param range = Tuple range of the spot prices we wish to examine
    # @param method = either Binomial Method 'bi' or Monte Carlo method 'mc'
    def __init__(self, S, r, sigma, T=1, range=(-10,1,10), method='bi'):
        
        ##Spot Price of asset
        self.S = S
        ##Risk free interest rate
        self.r = r
        ##Volatility of the asset (\f$\sigma\f$)
        self.sigma = sigma
        ##Time to expiry
        self.T = T
        ##Tuple range of the spot prices we wish to examine
        self.range = range
        ##either Binomial Method 'bi' or Monte Carlo method 'mc'
        self.meth = method
    
    ## Generates the option chain
    def generate(self):
        """Generates the option chain"""
        
        ##Dictionary of put values
        self.put = {}
        ##Dictionary of call values
        self.call = {}
        
        bottom = self.range[0]
        top = self.range[2]
        iter = self.range[1]
        
        if self.meth == 'bi':
            for E in range(bottom,top+iter,iter):
                solver = bi.binomial_euro(S=self.S,E=E,r=self.r,M=400,sigma=self.sigma,method='higham',T=self.T,opt='put')
                self.put[E] = solver.solve()
                del solver
            for E in range(bottom,top+iter,iter):
                solver = bi.binomial_euro(S=self.S,E=E,r=self.r,M=400,sigma=self.sigma,method='higham',T=self.T,opt='call')
                self.call[E] = solver.solve()
                del solver
        elif self.meth == 'mc':
            for E in range(bottom,top+iter,iter):
                solver = mc.mcfast_euro(S=self.S,E=E,r=self.r,M=20000,sigma=self.sigma,T=self.T,opt='put')
                self.put[E] = solver.solve()
                del solver
            for E in range(bottom,top+iter,iter):
                solver = mc.mcfast_euro(S=self.S,E=E,r=self.r,M=20000,sigma=self.sigma,T=self.T,opt='call')
                self.call[E] = solver.solve()
                del solver
        elif self.meth == 'all':
            pass
    
    ## Exports the result into a CSV file that may be opened in excel
    # @param filename = string filename to output. Please add extension CSV for proper functioning
    def export_csv(self,filename):
        f = open(filename,'w')
        
        f.write('Call, Days, Strike, Put\n')
        for key in sorted(self.call):
            f.write(str(self.call[key])+','+str(int(self.T*365))+','+str(key)+','+str(self.put[key])+'\n')
            
        
if __name__ == '__main__':
    #chain = option_chain(528.74,0.0011,0.3357,T=2.0/252.0,range=(480,5,580),method='bi')
    #chain.generate()
    #chain.export_csv('test.csv')