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

class option_chain():
    """Creates an option chain"""
    def __init__(self, S, r, sigma, T=1, range=(-10,1,10), method='bi'):
        self.S = S
        self.r = r
        self.sigma = sigma
        self.T = T
        self.range = range
        self.meth = method
    
    
    def generate(self):
        """Generates the option chain"""
        self.put = {}
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
    
    
    def export_csv(self,filename):
        f = open(filename,'w')
        
        f.write('Call, Days, Strike, Put\n')
        for key in sorted(self.call):
            f.write(str(self.call[key])+','+str(int(self.T*365))+','+str(key)+','+str(self.put[key])+'\n')
            
        
if __name__ == '__main__':
    chain = option_chain(528.74,0.0011,0.3357,T=2.0/252.0,range=(480,5,580),method='bi')
    chain.generate()
    chain.export_csv('test.csv')