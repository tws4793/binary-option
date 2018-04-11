import random
from math import exp, log, sqrt

import numpy as np
from numpy.linalg import inv

from scipy.stats import norm
from option import *

class BinaryOption(OptionPricing):
    """"A class to provide methods for binary option

    Attributes:
        S (float): Initial stock or index level
        K (float): Strike price
        T (float): Maturity (in year fraction)
        r (float): Constant risk-free short rate
        sigma (float): Volatility factor in diffusion term
    """
    bin_option_types = ['cash','asset']
    bin_asset_type = 1

    exp_cash = 0
    exp_asset = 0

    def __init__(self,S=0,K=0,r=0,sigma=0,T=0,q=0,asset_type=1):
        super().__init__(S,K,r,sigma,T,q)
        self.set_bin_asset_type(asset_type)

    # def __init__(self,S=0,K=0,r=0,sigma=0,T=0,q=0,asset_type=1):
    #     self.S = S
    #     self.K = K
    #     self.r = r
    #     self.sigma = sigma
    #     self.T = T
    #     self.q = q

    #     # Set the asset type at initialisation
    #     self.set_bin_asset_type(asset_type)
    
    # def set_variables(self,values):
    #     self.S = values[0]
    #     self.K = values[1]
    #     self.r = values[2]
    #     self.sigma = values[3]
    #     self.T = values[4]
    #     self.q = values[5]
    
    def set_bin_asset_type(self,set_type):
        """
        Attributes:
            set_type (string)
        """
        # Force variable into string
        set_type = str(set_type)

        # Checks if the value is in the array of option types
        if set_type not in self.bin_option_types:
            # Checks in the index is in the array of option types
            try:
                self.bin_option_types[int(set_type)-1]
                self.bin_asset_type = int(set_type)
            except:
                raise TypeError('Sorry, input for asset_type must be \'cash\' (1) or \'asset\' (2).')
        else:
            self.bin_asset_type = self.bin_option_types.index(set_type) + 1
        
        # Recalculate during each change
        self.exp_cash = self.get_cash_exp()
        self.exp_asset = self.get_asset_exp()
    
    def get_cash_exp(self):
        return exp(-self.r * self.T)
    
    def get_asset_exp(self):
        return exp(-self.q * self.T) * self.S

    def binary_black_scholes(self):
        d1 = (log(self.S / self.K) + (self.r - self.q + self.sigma ** 2 / 2) * self.T) / \
            (self.sigma * sqrt(self.T))
        d2 = (log(self.S / self.K) + (self.r - self.q - self.sigma ** 2 / 2) * self.T) / \
            (self.sigma * sqrt(self.T))
        
        # Calculate c and p
        factor = self.exp_cash if self.bin_asset_type == 1 else self.exp_asset
        d = d2 if self.bin_asset_type == 1 else d1
        
        c = factor * norm.cdf(d)
        p = factor * norm.cdf(-d)
        
        return (c,p)
    
    def binary_monte_carlo(self, simulations):
        def gbm():
            return self.S * np.exp((self.r - self.q - 0.5 * self.sigma**2) \
                * self.T + self.sigma \
                * np.sqrt(self.T) * random.gauss(0,1.0))
        
        def binary_payoff(mode, K, S_T):
            def call(K, S_T):
                return 1.0 if S_T >= K else 0.0
            def put(K, S_T):
                return 1.0 if S_T <= K else 0.0
            
            # The mode that was entered was not 1 (call) or 2 (put)
            if mode not in range(1,3):
                raise TypeError('Sorry, mode must be 1 (call) or 2 (put).')
            
            return call(K, S_T) if mode == 1 else put(K, S_T)
        
        payoffs_call = 0.0
        payoffs_put = 0.0

        for i in range(simulations):
            S_T = gbm()

            payoffs_call += binary_payoff(1, self.K, S_T)
            payoffs_put += binary_payoff(2, self.K, S_T)
        
        payoffs_avg = [payoffs_call / float(simulations),\
            payoffs_put / float(simulations)]
        factor = self.exp_cash if self.bin_asset_type == 1 else self.exp_asset

        c = factor * payoffs_avg[0]
        p = factor * payoffs_avg[1]

        return (c,p)