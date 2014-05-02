from pandas.io.data import DataReader
from pandas import np 
from numpy import *

## Calculate annualized historical volatility for a stock, defaulted to 10 days
#  @param sym Stock symbol for which we wish to get the historical volatility
#  @param days Time period in days we wish to calculate the volatility 
#
#  The calculation of annualized historical volatility is as follows: we begin by calculating the log returns of asset in question,
# \f[
#    U_i = \log \left ( \frac{ S(t_i) }{ S(t_{i-1}) } \right ),
#\f]
#over a predetermined period of time. It can be seen that this is the logarithm of the ratio of current and previous prices. In our context this is taken on a daily basis. We then calculate the mean of the returns
#\f[
#    a_M = \frac{ 1 }{ M } \sum_{i+1}^{M}(U_{n+1-i}).
#\f]
#Once the mean is calculated, we may then calculate the variance,
#\f[
#    b_M^2 = \frac{ 1 }{ M-1 }\sum_{i=1}^M (U_{n+1-i}-a_M)^2.
#\f]
#The general form of annualized volatility is stated to be,
#\f[
#    \sigma_a = \sqrt{Tb^2_M},
#\f]
#Where \f$T\f$ is a timeframe in sample space, which we early stated was daily for our purposes. Since we are considering annualized volatility for a year time-frame, it is known there are 252 trading days in a year and this equation becomes
#\f[
#    \sigma_a = \sqrt{252b^2_M}.
#\f]
def hist_vol(sym, days=10):
    try:
        quotes = DataReader(sym, 'yahoo')['Close'][-days:]
    except Exception:
        print "Problem getting historical volatility!"
        raise SystemExit(code)
        return None, None
    logreturns = np.log(quotes / quotes.shift(1))
    vol = np.sqrt(252*logreturns.var()) #252 trading days in year (annualized volatility)
    return float(vol)

    

if __name__ == '__main__':
    #print hist_vol('goog',20)