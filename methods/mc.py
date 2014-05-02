## @package pyFi.methods.mc
# Contains Monte Carlo method solvers for option valuation

from numpy import *
from numpy.random import standard_normal as rm
from numpy.random import randn


## Monte Carlo Method solver for European vanilla options
class montecarlo_euro():
    def __init__(self, S, E, r, sigma, T=1, M=10000, opt='call'):
        ## Called upon initialization of the Monte Carlo method for Europeans
        #
        # @param S = Spot Price of asset
        # @param E = Exercise / Strike price
        # @param r = Risk free interest rate
        # @param sigma = Volatility of the asset (\f$\sigma\f$)
        # @param T = Time to expiry (defaults to 1.0)
        # @param M = Number of random samples
        # @param opt = 'call' or 'put' (defaults to 'put')
        self.S = S
        self.E = E
        self.sigma = sigma
        self.r = r
        self.T = T
        self.M = M
        self.opt = opt
    
    ## Produces the Monte Carlo solution 
    #
    # \f[e^{-rT}\Lambda \left ( S_0 \exp \left [ \left ( r - \frac{ 1 }{ 2 }\sigma^2 \right )T + \sigma \sqrt{T}Z \right ] \right )\f]
    # where
    # \f[Z \sim N(0,1)\f]
    def solve(self):
        """comment"""
        V = zeros(self.M)
        
        for i in range(1,self.M):
            Sf = self.S*exp((self.r-0.5*self.sigma*self.sigma)*self.T+self.sigma*sqrt(self.T)*rm())
            V[i] = exp(-self.r*self.T)*self.gamma(Sf)
        self.meany = mean(V)
        self.bet = std(V)
        return self.meany
    
    ## Evaluates whether option or put
    def gamma(self,sf):
        if self.opt=='put':
            gam = max(self.E-sf,0)
        elif self.opt=='call':
            gam = max(sf-self.E,0)
        return gam
    
    ## Returns the certainty bars of the solution
    def certainty(self):
        return (self.meany-self.bet,self.meany+self.bet)
    
## Vectorized fast Monte Carlo Method solver for European vanilla options
class mcfast_euro():
    ## Called upon initialization of the Monte Carlo method for Europeans
    #
    # @param S = Spot Price of asset
    # @param E = Exercise / Strike price
    # @param r = Risk free interest rate
    # @param sigma = Volatility of the asset (\f$\sigma\f$)
    # @param T = Time to expiry (defaults to 1.0)
    # @param M = Number of random samples
    # @param opt = 'call' or 'put' (defaults to 'put')
    def __init__(self, S, E, r, sigma, T=1, M=10000, opt='call'):
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
    
    ## Vectorized Solution
    #
    # We initialize a random array of variables \f$Z\f$ then calculate the constants
    # \f[ C_1 = \left ( r - \frac{1}{2} \sigma^2 \right )T\f]
    # \f[ C_2 = \sigma \sqrt{T} \f]
    # \f[ C_3 = e^{-rt} \f]
    # We then calculate
    # \f[\boldsymbol S_f = S \circ \exp \left ( C_1 + C_2 \circ \boldsymbol Z \right )\f]
    # and yield the answer
    # \f[ W(S) = \overline{ C_3 \circ \boldsymbol \Lambda \left ( S_f \right ) }\f]
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