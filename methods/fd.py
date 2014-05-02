## @package pyFi.methods.fd
# Contains finite difference approaches for the Black-Scholes PDE

from numpy import *
from numpy.linalg import solve

## Implicit method for solving BS Equation
# @todo: DOES NOT CURRENTLY WORK
# @todo: REQUIRES DOCUMENTATION
class BS_fd_implicit():
    ## Called upon initialization of the Finite Difference method for Europeans
    #
    # @param S = Spot Price of asset
    # @param E = Exercise / Strike price
    # @param r = Risk free interest rate
    # @param sigma = Volatility of the asset (\f$\sigma\f$)
    # @param L = Length of spot domain
    # @param T = Time to expiry (defaults to 1.0)
    # @param k = Number of timesteps
    # @param M = Discretization points
    # @param opt = 'call' or 'put' (defaults to 'put')
    def __init__(self, S, E, r, sigma, L, T=1.0, k=400, M=400, opt='put'):
        """Initializes all class variables and calculates constants"""
        self.S = S
        self.E = E
        self.r = r
        self.sigma = sigma
        self.L = L
        self.T = T
        self.k = k
        self.M = M
        self.opt = opt
        
        self.dt = T/k
        self.h  = float(L)/(float(M)-1.0)
        
        self.n = array([x for x in range(0,self.M)])
        
        self.Mat = self.build_LHS_matrix()
    
    
    def build_LHS_matrix(self):
        """Builds the LHS operator matrix (A) for the Black Scholes equation of the form Ax=b"""
        a = 0.5*self.dt*(self.sigma**2 * self.n**2 - self.r*self.n)
        b = 1.0 - self.dt*(self.sigma**2 * self.n**2 + self.r)
        c = 0.5*self.dt*(self.sigma**2*self.n**2 + self.r*self.n)
        
        Mat = diag(b) + diag(c[0:self.M-1],1) + diag(a[1:self.M],-1)
        Mat[0,:] = 0
        Mat[0,0] = 1
        Mat[self.M-1,:] = 0
        Mat[self.M-1,self.M-1] = 1
        return Mat
    
    
    def strike(self):
        """creates a spot array"""
        return self.h*self.n 
    
    
    def init_values_call(self):
        """Initial condition for a call option"""
        val = maximum(self.strike()-self.E,0)
        return reshape(val,(self.M,1))
    
    
    def init_values_put(self):
        """Initial condition for a put option"""
        val = maximum(self.E-self.strike(),0)
        return reshape(val,(self.M,1))
    
    
    def update_bc(self,U,i):
        """Updates the explicit boundary condition for a call or a put"""
        if self.opt=='put':
            U[0,0] = self.E*exp(-self.r*self.dt*(i+1))
            U[self.M-1,0] = 0
        elif self.opt=='call':
            U[self.M-1,0] = self.L-self.E*exp(-self.r*self.dt*(i+1))
            U[0,0] = 0
        return U
    
    
    def solve_1d_surface(self):
        """Returns a 1D space line for the final time of BS equation"""
        if self.opt=='call':
            RHS = self.init_values_call()
        elif self.opt=='put':
            RHS = self.init_values_put()
        
        for i in range(self.k):
            U = solve(self.Mat,RHS)
            U = self.update_bc(U,i)
            RHS = U
        
        return U
    
    
    def solve_2d_surface(self):
        """Returns a time-space surface of the solution of the BS equation"""
        sol = []
        
        if self.opt=='call':
            Uold = self.init_values_call()
        elif self.opt=='put':
            Uold = self.init_values_put()
        
        for i in range(self.k):
            U = dot(self.Mat,Uold)
            U = self.update_bc(U,i)
            Uold = U
            sol.append(U)
        
        return sol
    
    
    def interp_solution(self,u,i0,i1):
        """Interpolates the solution from the 1d surface for the strike price we are interested in"""
        x = self.n*self.h
        x0 = x[i0]
        x1 = x[i1]
        y1 = u[i1]
        y0 = u[i0]
        
        y = y0 + (y1-y0)*(self.S-float(x0))/(float(x1)-float(x0))
        return y
    
    
    def solve(self):
        """Solves for a particular value by calling interpolation routine"""
        U = self.solve_1d_surface()
        where = self.S / self.h
        top = int(ceil(where))
        bottom = int(floor(where))
        
        if top==bottom:
            return U[top]
        else:
            return self.interp_solution(U,bottom,top)

## Solution for the Black Scholes equation using an explicit central-space finite difference solver
#
# Solves the Black-Scholes equation,
# \f[\boldsymbol V_t + \frac{1 }{ 2 } \sigma^2 \boldsymbol S^2 \boldsymbol V_{SS} + r \boldsymbol S \boldsymbol V_S - r \boldsymbol V = 0\f]
# with boundary conditions for a put
# \f[V(0,\tau) = Ee^{-r\tau}\f]
# \f[V(L,\tau) = 0\f]
# or for a call
# \f[V(0,\tau) = 0\f]
# \f[V(L,\tau) = S-Ee^{-r\tau}\f]
class BS_fd_explicit():
    ## Called upon initialization of the Finite Difference method for Europeans
    #
    # @param S = Spot Price of asset
    # @param E = Exercise / Strike price
    # @param r = Risk free interest rate
    # @param sigma = Volatility of the asset (\f$\sigma\f$)
    # @param L = Length of spot domain
    # @param T = Time to expiry (defaults to 1.0)
    # @param k = Number of timesteps
    # @param M = Discretization points
    # @param opt = 'call' or 'put' (defaults to 'put')
    def __init__(self, S, E, r, sigma, L, T=1.0, k=400, M=400, opt='put'):
        self.S = S
        self.E = E
        self.r = r
        self.sigma = sigma
        self.L = L
        self.T = T
        self.k = k
        self.M = M
        self.opt = opt
        
        self.dt = T/k
        self.h  = float(L)/(float(M)-1.0)
        
        self.n = array([x for x in range(0,self.M)])
        
        self.Mat = self.build_matrix()
        
    
    ## Builds the operator matrix for explicit FD
    #
    # \f[ \begin{bmatrix}1 & 0 & 0 & \cdots & 0 & 0 & 0 \\ a & b & c & \cdots & 0 & 0 & 0 \\ 0 & a & b &        & 0 & 0 & 0 \\ \vdots  & \vdots &   & \ddots  & & \vdots & \vdots  \\ 0 & 0 & 0 & \cdots & b & c & 0 \\ 0 & 0 & 0 & \cdots & a & b & c \\ 0 & 0 & 0 & \cdots & 0 & 0 & 1 \end{bmatrix}\f]
    #
    # Where:
    # 
    # \f[ a = \frac{ 1 }{ 2 }(\Delta t)rn[\sigma^2 n - r]\\
    #  b = 1- (\Delta t)rn[\sigma^2 n^2 - r]\\
    #  = \frac{ 1 }{ 2 }(\Delta t)rn [ \sigma^2 n + r]. \f]
    def build_matrix(self):
        a = 0.5*self.dt*(self.sigma**2 * self.n**2 - self.r*self.n)
        b = 1.0 - self.dt*(self.sigma**2 * self.n**2 + self.r)
        c = 0.5*self.dt*(self.sigma**2*self.n**2 + self.r*self.n)
        
        Mat = diag(b) + diag(c[0:self.M-1],1) + diag(a[1:self.M],-1)
        Mat[0,:] = 0
        Mat[0,0] = 1
        Mat[self.M-1,:] = 0
        Mat[self.M-1,self.M-1] = 1
        return Mat
    
    ## Creates the spatial array
    def strike(self):
        return self.h*self.n 
    
    ## Initial condition for a call option
    #
    # \f[V(S,\tau) = \max(E-S,0) \f]
    def init_values_call(self):
        val = maximum(self.strike()-self.E,0)
        return reshape(val,(self.M,1))
    
    ## Initial condition for a put option
    #
    # \f[ V(S,\tau) = \max(E-S,0) \f]
    def init_values_put(self):
        val = maximum(self.E-self.strike(),0)
        return reshape(val,(self.M,1))
    
    ## Updates the boundary condition for each timestep
    # 
    # for a put
    # \f[V(0,\tau) = Ee^{-r\tau}\f]
    # \f[V(L,\tau) = 0\f]
    # or for a call
    # \f[V(0,\tau) = 0\f]
    # \f[V(L,\tau) = S-Ee^{-r\tau}\f]
    def update_bc(self,U,i):
        if self.opt=='put':
            U[0,0] = self.E*exp(-self.r*self.dt*(i+1))
        elif self.opt=='call':
            U[self.M-1,0] = self.L-self.E*exp(-self.r*self.dt*(i+1))
        return U
    
    ## Solves the equation and returns the entire 1d surface
    def solve_1d_surface(self):
        if self.opt=='call':
            Uold =  self.init_values_call()
        elif self.opt=='put':
            Uold = self.init_values_put()
        
        for i in range(self.k):
            U = dot(self.Mat,Uold)
            U = self.update_bc(U,i)
            Uold = U
        
        return U
    
    ## Solves the equation and returns a surface as a function of time and space
    def solve_2d_surface(self):
        """Returns a time-space surface of the solution of the BS equation"""
        sol = []
        
        if self.opt=='call':
            Uold = self.init_values_call()
        elif self.opt=='put':
            Uold = self.init_values_put()
        
        for i in range(self.k):
            U = dot(self.Mat,Uold)
            U = self.update_bc(U,i)
            Uold = U
            sol.append(U)
        
        return sol
    
    ## Linear Interpolant routine
    #
    # \f[y = y_0 + (y_1-y_0)\frac{ x-x_0 }{ x_1-x_0 }\f]
    def interp_solution(self,u,i0,i1):
        """Interpolates the solution from the 1d surface for the strike price we are interested in"""
        x = self.n*self.h
        x0 = x[i0]
        x1 = x[i1]
        y1 = u[i1]
        y0 = u[i0]
        
        y = y0 + (y1-y0)*(self.S-float(x0))/(float(x1)-float(x0))
        return y
    
    ## This method solves the eactual problem for a single spot value we are intestested in
    #
    # Solves the 1-d surface for the final time then calls the interpolation routine to determine
    # the value at the spot price of interest
    def solve(self):
        """Solves for a particular value by calling interpolation routine"""
        U = self.solve_1d_surface()
        where = self.S / self.h
        top = int(ceil(where))
        bottom = int(floor(where))
        
        if top==bottom:
            return U[top]
        else:
            return self.interp_solution(U,bottom,top)
    
    


if __name__ == '__main__':
    solver = BS_fd_explicit(S=4,E=5,r=0.04,sigma=0.3,M=20,L=10,k=25000,opt='put')
    #print solver.solve_1d_surface()
    print solver.solve()
    
    #solver2 = BS_fd_implicit(S=3,E=2,r=0.05,sigma=0.3,M=6,L=10,k=100,opt='put')
    #print solver2.solve_1d_surface()
    #print solver2.solve()