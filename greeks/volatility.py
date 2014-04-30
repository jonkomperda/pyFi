from pandas.io.data import DataReader
from pandas import np 
from numpy import *

 
def hist_vol(sym, days=10):
    "Calculate annualized historical volatility for a stock, defaulted to 10 days"
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
    print hist_vol('goog',20)