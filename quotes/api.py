import urllib
import re

## Gets quote for a symbol from Google Finance
# @param symbol= String name of stock symbol. For example 'GOOG' would be Google
#
# Uses urllib to get quote from http://finance.google.com/finance?q=
def get_quote_google(symbol):
    base_url = 'http://finance.google.com/finance?q='
    content = urllib.urlopen(base_url + symbol).read()
    m = re.search('id="ref_.*?_l".*?>(.*?)<', content)
    if m:
        quote = m.group(1)
    else:
        print 'No quote available for: ' + symbol
        quote = float(raw_input('Please enter symbol quote: $'))
    return float(quote)

## Gets quote for a symbol from Yahoo Finance
# @param symbol= String name of stock symbol. For example 'GOOG' would be Google
#
# Uses urllib to get quote from http://finance.yahoo.com/q?s=
def get_quote_yahoo(symbol):
    base_url = 'http://finance.yahoo.com/q?s='
    content = urllib.urlopen(base_url + symbol).read()
    m = re.search('id="yfs_l84_.*?>(.*?)<',content)
    if m:
        quote = m.group(1)
    else:
        print 'No quote available for: ' + symbol
        quote = float(raw_input('Please enter symbol quote: $'))
    return float(quote)

## Gets the risk free rate
# Uses 1 year Treasury Note rate from bankrate.com
def get_risk_free():
    base_url = 'http://www.bankrate.com/rates/interest-rates/1-year-treasury-rate.aspx'
    content = urllib.urlopen(base_url).read()
    m = re.search('class="tabledataoddnew">0(.*?)<',content)
    if m:
        quote = float(m.group(1))/100.0
    else:
        print 'No T-Note data available...'
        quote = float(raw_input('Please enter risk free interest rate: '))
    return quote
    

if __name__ == '__main__':
    #print get_quote_google('GOOG')
    #print get_quote_yahoo('AAPL')
    print get_risk_free()