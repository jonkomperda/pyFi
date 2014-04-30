from numpy import *
from numpy.random import standard_normal as rm
from numpy.random import randn

class montecarlo_euro():
    """docstring for montecarlo"""
    def __init__(self, S, E, r, sigma, T=1, M=10000, opt='call'):
        self.S = S
        self.E = E
        self.sigma = sigma
        self.r = r
        self.T = T
        self.M = M
        self.opt = opt
    
    
    def solve(self):
        """comment"""
        V = zeros(self.M)
        
        for i in range(1,self.M):
            Sf = self.S*exp((self.r-0.5*self.sigma*self.sigma)*self.T+self.sigma*sqrt(self.T)*rm())
            V[i] = exp(-self.r*self.T)*self.gamma(Sf)
        self.meany = mean(V)
        self.bet = std(V)
        return self.meany
    
    
    def gamma(self,sf):
        if self.opt=='put':
            gam = max(self.E-sf,0)
        elif self.opt=='call':
            gam = max(sf-self.E,0)
        return gam
    
    
    def certainty(self):
        return (self.meany-self.bet,self.meany+self.bet)
    

class mcfast_euro():
    """docstring for mcfast"""
    def __init__(self, S, E, r, sigma, dt=0.0001, T=1, M=10000, opt='call'):
        self.S = S
        self.E = E
        self.sigma = sigma
        self.r = r
        self.T = T
        self.M = M
        self.opt = opt
    
    
    def gamma(self,sf):
        """docstring for gamma"""
        if self.opt=='put':
            gam = maximum(self.E-sf,0)
        if self.opt=='call':
            gam = maximum(sf-self.E,0)
        return gam
    
    
    def solve(self):
        randarray = randn(self.M)
        
        c1 = (self.r-0.5*self.sigma*self.sigma)*self.T
        c2 = self.sigma*sqrt(self.T)
        
        Sf = self.S*exp(c1+c2*randarray)
        
        c3 = exp(-self.r*self.T)
        
        gam = self.gamma(Sf)
        
        V  = c3*gam
        
        meany = mean(V)
        
        return meany
        
        
    
if __name__ == '__main__':
    #solve = montecarlo_euro(4,5,0.04,0.3,M=1000000,opt='put')
    #print solve.solve()
    
    solve2 = mcfast_euro(4,5,0.04,0.3,M=500000000,opt='put')
    print solve2.solve()