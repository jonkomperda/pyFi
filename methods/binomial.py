from numpy import *
class binomial_euro():
    def __init__(self, S, E, r, sigma, T=1.0, M=400, p=0.5, opt='call', method='higham'):
        self.S = S
        self.E = E
        self.r = r
        self.T = T
        self.sigma = sigma
        self.M = M
        self.p = p
        self.dt = T/M
        
        if method=='higham':
            self.u,self.d = self.ud_higham()
        elif method=='crr':
            self.u,self.d = self.ud_crr()
        elif method=='kwok':
            self.u,self.d = self.ud_kwok()
        else:
            self.u,self.d = 0
        
        if opt=='call':
            self.W = self.final_values_call()
        elif opt=='put':
            self.W = self.final_values_put()
    
    
    def ud_higham(self):
        """docstring for ud_higham"""
        u = exp(self.sigma*sqrt(self.dt) + (self.r-0.5*self.sigma**2.0)*self.dt)
        d = exp(-self.sigma*sqrt(self.dt) + (self.r-0.5*self.sigma**2.0)*self.dt) 
        return u,d
    
    
    def ud_crr(self):
        """ comment """
        u = exp(self.sigma*sqrt(self.dt))
        d = exp(-self.sigma*sqrt(self.dt))
        return u,d
    
    
    def ud_kwok(self):
        """docstring for ud_kwok"""
        u = exp(self.sigma*sqrt(self.dt)*(1.0+sqrt(exp(self.sigma*self.sigma*self.dt-1.0))))
        d = exp(self.sigma*sqrt(self.dt)*(1.0-sqrt(exp(self.sigma*self.sigma*self.dt-1.0))))
        return u,d
    
    
    def tree(self):
        """docstring for tree"""
        m1 = range(self.M,-1,-1)
        m2 = range(0,self.M+1)
        dp = self.d**m1
        up = self.u**m2
        return dp,up
    
    
    def final_values_call(self):
        """docstring for final_values"""
        dp,up=self.tree()
        #w = zeros(self.M+1)
        w = maximum(self.S*dp*up-self.E,0)
        return w
    
    
    def final_values_put(self):
        """docstring for final_values_put"""
        dp,up=self.tree()
        #w = zeros(self.M+1)
        w = maximum(self.E-self.S*dp*up,0)
        return w
    
    
    def solve(self):
        """docstring for solve"""
        w = copy(self.W)
        for k in range(self.M,0,-1):
            w = exp(-self.r*self.dt)*(self.p*w[1:k+1] + (1.0-self.p)*w[0:k])
        return float(w)
    
    

class binomial_amer(binomial_euro):
    """docstring for binomial_amer"""
    def __init__(self, S,E,r,sigma,T=1.0,M=400,p=0.5,opt='call',method='higham'):
        self.S=S
        self.E=E
        self.r=r
        self.sigma=sigma
        self.T=T
        self.M=M
        self.p=p
        self.dt = T/M
        
        if method=='higham':
            self.u,self.d = self.ud_higham()
        elif method=='crr':
            self.u,self.d = self.ud_crr()
        elif method=='kwok':
            self.u,self.d = self.ud_kwok()
        else:
            self.u,self.d=0
        
        self.dp,self.up = self.tree()
        
        if opt=='call':
            self.W = self.final_values_call()
        elif opt=='put':
            self.W = self.final_values_put()
    
    
    def solve(self):
        """docstring for solve"""
        w = copy(self.W)
        for k in range(self.M,0,-1):
            temp = self.S*self.dp[self.M-i+1:self.M]*self.up[0:i]
            w = maximum(maximum(self.E-temp,0),exp(-self.r*self.dt)*(self.p*w[1:i]+(1-self.p)*w[0:i]))
        return float(w)
    


if __name__ == '__main__':
    solve = binomial_euro(3,2,0.05,0.3,M=40000,method='higham')
    #print solve.u
    #print solve.d
    #print solve.W
    print solve.solve()