import numpy as np
import pandas as pd
import statistics
from config import *

def get_z(N):
    z = np.random.standard_normal(N)
    z = np.append(z,-z)
    return z

def get_ticker_data(ticker):
    def get_sig(data):
        log = []
        rng = 63
        x = data[:rng].iloc[0]['Adj. Close']
        for i in range(1, rng):
            y = data[:rng].iloc[i]['Adj. Close']
            log.append(np.log(y/x))
            x = y
        sigma = statistics.stdev(log)*np.sqrt(252)
        return sigma
    
    ext = 'csv'
    
    url_quandl = 'https://www.quandl.com/api/v3/datasets/WIKI/' + ticker + '.' + ext + '?api-key=' + QUANDL_API_KEY
    url_yahoo = 'https://finance.yahoo.com/quote/' + ticker + '/key-statistics'
    url_treasury = 'https://'

    data_quandl = pd.read_csv(url_quandl, index_col=0)
    data_yahoo = pd.read_html(url_yahoo, index_col=0)
    # data_treasury = pd.read_(url_treasury, index_col=0)

    return {
        'S': float(data_quandl['Adj. Close'][0]), # Latest Adjusted closing
        'r': 0.017, # 3 Month US T-Bill
        'sigma': get_sig(data_quandl), # Historical Stock Data
        'q': float(data_yahoo[9][1][3].strip('%'))/100, # Trailing Annual Dividend Yield
    }