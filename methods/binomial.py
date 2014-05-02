## @package pyFi.methods.binomial
# Contains binomial method solvers for options
from numpy import *
class binomial_euro():
    ## Called upon initialization of the binomial method for Europeans
    #
    # @param S = Spot Price of asset
    # @param E = Exercise / Strike price
    # @param r = Risk free interest rate
    # @param sigma = Volatility of the asset (\f$\sigma\f$)
    # @param T = Time to expiry (defaults to 1.0)
    # @param M = Size of binomial tree (defaults to 400)
    # @param p = probability (defaults to 0.4)
    # @param opt = 'call' or 'put' (defaults to 'call')
    # @param method = u,d method, either 'higham', 'kwok' or 'crr' (defaults to 'higham')
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
    
    ## u,d values from Higham
    #
    # \f[ u = e^{\sigma \sqrt{\delta t} + \left ( r - \frac{ \sigma^2 }{ 2 } \right ) \delta t}\f]
    # \f[ d = e^{-\sigma \sqrt{\delta t} + \left ( r - \frac{ \sigma^2 }{ 2 } \right ) \delta t} \f]
    def ud_higham(self):
        u = exp(self.sigma*sqrt(self.dt) + (self.r-0.5*self.sigma**2.0)*self.dt)
        d = exp(-self.sigma*sqrt(self.dt) + (self.r-0.5*self.sigma**2.0)*self.dt) 
        return u,d
    
    ## u,d from Cox, Ross, and Rubenstein
    #
    # \f[ u  = e^{\sigma \sqrt{\delta t}} \f]
    # \f[ d = e^{\sigma \sqrt{-\delta t}} \f]
    def ud_crr(self):
        u = exp(self.sigma*sqrt(self.dt))
        d = exp(-self.sigma*sqrt(self.dt))
        return u,d
    
    ## u,d values from Kwok
    #
    # \f[ u = e^{\sigma \sqrt{\delta t}} \left ( 1 + \sqrt{e^{\sigma^2 \delta t} - 1} \right )\f]
    # \f[ d = e^{\sigma \sqrt{\delta t}} \left ( 1 - \sqrt{e^{\sigma^2 \delta t} - 1} \right )\f]
    def ud_kwok(self):
        u = exp(self.sigma*sqrt(self.dt)*(1.0+sqrt(exp(self.sigma*self.sigma*self.dt-1.0))))
        d = exp(self.sigma*sqrt(self.dt)*(1.0-sqrt(exp(self.sigma*self.sigma*self.dt-1.0))))
        return u,d
    
    ## Builds the tree
    #
    # \f[M_1 = [\,M \quad M-1 \quad M-2 \quad ... \quad 0 \, ] \f]
    #\f[M_2  = [\,0 \quad 1 \quad 2 \quad ... \quad M \, ],\f]
    #\f[d_p  = d^{M_1}\f]
    #\f[u_p  = u^{M_2}\f]
    def tree(self):
        m1 = range(self.M,-1,-1)
        m2 = range(0,self.M+1)
        dp = self.d**m1
        up = self.u**m2
        return dp,up
    
    ## Determines the final values for a call
    #
    # \f[ \boldsymbol W = \max(S \circ d_p \circ u_p - E, 0) \quad \text{for a call} \f]
    def final_values_call(self):
        """docstring for final_values"""
        dp,up=self.tree()
        #w = zeros(self.M+1)
        w = maximum(self.S*dp*up-self.E,0)
        return w
    
    ## Determines the final values for a put
    #
    # \f[ \boldsymbol W = \max(E-S \circ d_p \circ u_p, 0), \quad \text{for a put} \f]
    def final_values_put(self):
        dp,up=self.tree()
        #w = zeros(self.M+1)
        w = maximum(self.E-self.S*dp*up,0)
        return w
    
    ## Calculates the recursion for solving a binomial tree
    #
    # \f[ V_n^i = e^{-r\delta t}\left ( pV_{n+1}^{i+1} + (1-p) V_n^{i+1} \right ), quad 0 \le n \le i, \quad 0 \le i \le M-1 \f]
    # Returns W, the asset price.
    def solve(self):
        w = copy(self.W)
        for k in range(self.M,0,-1):
            w = exp(-self.r*self.dt)*(self.p*w[1:k+1] + (1.0-self.p)*w[0:k])
        return float(w)
    
    

class binomial_amer(binomial_euro):
    ## Called upon initialization of the binomial method for Europeans
    #
    # @param S = Spot Price of asset
    # @param E = Exercise / Strike price
    # @param r = Risk free interest rate
    # @param sigma = Volatility of the asset (\f$\sigma\f$)
    # @param T = Time to expiry (defaults to 1.0)
    # @param M = Size of binomial tree (defaults to 400)
    # @param p = probability (defaults to 0.4)
    # @param opt = 'call' or 'put' (defaults to 'call')
    # @param method = u,d method, either 'higham', 'kwok' or 'crr' (defaults to 'higham')
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
    
    ## Solves the recurion relation for an American Option
    def solve(self):
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