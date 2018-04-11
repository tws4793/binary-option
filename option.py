import random
from math import exp, log, sqrt

import numpy as np
from numpy.linalg import inv

from scipy.stats import norm

class OptionPricing:
    """"A class to provide methods for option pricing

    Attributes:
        S (float): Initial stock or index level
        K (float): Strike price
        T (float): Maturity (in year fraction)
        r (float): Constant risk-free short rate
        sigma (float): Volatility factor in diffusion term
    """

    def __init__(self,S=0,K=0,r=0,sigma=0,T=0,q=0):
        self.S = S
        self.K = K
        self.r = r
        self.sigma = sigma
        self.T = T
        self.q = q
    
    def set_variables(self,values):
        self.S = values[0]
        self.K = values[1]
        self.r = values[2]
        self.sigma = values[3]
        self.T = values[4]
        self.q = values[5]

    def black_scholes(self):
        d1 = (log(self.S/self.K)+(self.r-self.q+self.sigma**2/2)*self.T)/(self.sigma*sqrt(self.T))
        d2 = d1-self.sigma*sqrt(self.T)
        
        c = self.S*exp(-self.q*self.T)*norm.cdf(d1)-self.K*exp(-self.r*self.T)*norm.cdf(d2)
        p = self.K*exp(-self.r*self.T)*norm.cdf(-d2)-self.S*exp(-self.q*self.T)*norm.cdf(-d1)
        
        return (c,p)

    def binomial_tree(self,N=100):
        # 1.
        dT=self.T/N
        u=np.exp(self.sigma*np.sqrt(dT))
        d=1/u
        p=(np.exp((self.r-self.q)*dT)-d)/(u-d)

        # 2. 
        fc=np.zeros((N+1,N+1))
        fp=np.zeros((N+1,N+1))
        j=np.arange(0,N+1,1)
        Ss=self.S*(u**j)*(d**(N-j))
        fc[N,:]=np.maximum(0, Ss-self.K)
        fp[N,:]=np.maximum(0, self.K-Ss)

        # 3.
        p1=1-p
        ert=np.exp(-self.r*dT)
        for i in range(N-1,0-1,-1):
            fc[i,0:i+1]=ert*(p*fc[i+1,0+1:i+1+1]+p1*fc[i+1,0:i+1]) 
            fp[i,0:i+1]=ert*(p*fp[i+1,0+1:i+1+1]+p1*fp[i+1,0:i+1])
        
        # 4.
        c=fc[0][0]
        p=fp[0][0]

        return (c, p)
    
    def monte_carlo(self,z):
        #2.
        ST=self.S*np.exp((self.r-self.q-0.5*self.sigma**2)*self.T+self.sigma*np.sqrt(self.T)*z)
        
        #3.
        hTc=np.maximum(ST-self.K,0)
        hTp=np.maximum(self.K-ST,0)

        #4.
        c=np.exp(-self.r*self.T)*np.mean(hTc) 
        p=np.exp(-self.r*self.T)*np.mean(hTp)

        return (c, p)
    
    def implicit(self,M=100,N=1000):
        #1.
        dt= self.T/N
        Smax= 2*self.K
        dS= Smax/M
        #2
        Fc=np.zeros((M+1,N+1)) #Call Option
        Fp=np.zeros((M+1,N+1)) #Put option
        for j in range(M+1):
            Fc[j][N]=np.maximum(j*dS-self.K,0)
            Fp[j][N]=np.maximum(self.K-j*dS,0)   

        A=np.zeros((M+1,M+1)) 
        A[M,M]=1
        A[0,0]=1
        for j in range(1,M):
            A[j][j-1]=0.5*dt*((self.r-self.q) * j - self.sigma**2 * j**2)
            A[j][j]= 1+dt*(self.sigma**2 * j**2 + self.r)
            A[j][j+1]= -0.5*dt*(self.sigma**2 * j**2 + (self.r-self.q) * j)
        #3.1, 3.2    
        for i in range(N-1,-1,-1):
            #Call Option
            FcN= Fc[:,[i+1]] #Get the Fi+1
            FcN[0][0] =0
            FcN[M][0]= Smax - (self.K* np.exp(-self.r * (N - i) *dt))
            Fc[:,[i]]= np.dot(inv(A),FcN)
            #PutOption
            FpN= Fp[:,[i+1]]
            FpN[0][0] = self.K* np.exp(-self.r * (N - i) *dt)
            FpN[M][0]= 0 
            Fp[:,[i]]= np.dot(inv(A),FpN)
        #4
        k= np.int32(np.floor(self.S/dS))
        #5
        C = Fc[k][0] + (Fc[k+1][0]-Fc[k][0])/dS * (self.S - k * dS)
        P = Fp[k][0] + (Fp[k+1][0]-Fp[k][0])/dS * (self.S - k * dS)
        return (C, P)

    def explicit(self,M=100,N=1000):
        #1.
        dt= self.T/N
        Smax= 2*self.K
        dS= Smax/M
        #2
        Fc=np.zeros((M+1,N+1))
        Fp=np.zeros((M+1,N+1))
        for j in range(M+1):
            Fc[j][N]=np.maximum(j*dS-self.K,0)
            Fp[j][N]=np.maximum(self.K-j*dS,0)   
        #3.1
        A=np.zeros((M+1,M+1)) 
        A[M,M]=1
        A[0,0]=1
        for j in range(1,M):
            A[j][j-1]= 0.5*dt*(self.sigma**2 * j**2 - (self.r-self.q) * j)
            A[j][j]= 1-dt*(self.sigma**2 * j**2 + self.r)
            A[j][j+1]= 0.5*dt*(self.sigma**2 * j**2 + (self.r-self.q) * j)
        #3.2
        for i in range(N-1,-1,-1):
            #Call Option
            FcN= Fc[:,[i+1]] #Get the Fi+1
            Fc[:,[i]]= np.dot(A,FcN) #Update function with matrix multiplication
            Fc[0][i] =0
            Fc[M][i]= Smax - (self.K* np.exp(-self.r * (N - i) *dt)) 
            #PutOption
            FpN= Fp[:,[i+1]]
            Fp[:,[i]]= np.dot(A,FpN)
            Fp[0][i] = self.K* np.exp(-self.r * (N - i) *dt)
            Fp[M][i]= 0  
        #4
        k= np.int32(np.floor(self.S/dS))
        #5
        C = Fc[k][0] + (Fc[k+1][0]-Fc[k][0])/dS * (self.S - k * dS)
        P = Fp[k][0] + (Fp[k+1][0]-Fp[k][0])/dS * (self.S - k * dS)
        #print(A)
        #print(Fc)
        #print(Fp)
        return (C, P)
    
    def crank_nicolson(self,M=100,N=1000):
        #1.
        dt= self.T/N
        Smax= 2*self.K
        dS= Smax/M
        #2
        Fc=np.zeros((M+1,N+1)) #Call Option
        Fp=np.zeros((M+1,N+1)) #Put option
        for j in range(M+1):
            Fc[j][N]=np.maximum(j*dS-self.K,0)
            Fp[j][N]=np.maximum(self.K-j*dS,0)   

        M1=np.zeros((M+1,M+1)) 
        M1[M,M]=1
        M1[0,0]=1
        for j in range(1,M):
            M1[j][j-1]= -0.25*dt*(self.sigma**2 * j**2 - (self.r-self.q) * j)
            M1[j][j]= 1-(-0.5*dt*(self.sigma**2 * j**2 + self.r))
            M1[j][j+1]= -0.25*dt*(self.sigma**2 * j**2 + (self.r-self.q) * j)
        M2=np.zeros((M+1,M+1)) 
        M2[M,M]=1
        M2[0,0]=1
        for j in range(1,M):
            M2[j][j-1]= 0.25*dt*(self.sigma**2 * j**2 - (self.r-self.q) * j)
            M2[j][j]= 1+(-0.5*dt*(self.sigma**2 * j**2 + self.r))
            M2[j][j+1]= 0.25*dt*(self.sigma**2 * j**2 + (self.r-self.q) * j)

        #3.1, 3.2    
        for i in range(N-1,-1,-1):
            #Call Option
            FcN= np.dot(M2, Fc[:,[i+1]]) #Get the Fi+1
            FcN[0][0] =0
            FcN[M][0]= Smax - (self.K* np.exp(-self.r * (N - i) *dt))
            Fc[:,[i]]= np.dot(inv(M1),FcN)
            #PutOption
            FpN= np.dot(M2, Fp[:,[i+1]])
            FpN[0][0] = self.K* np.exp(-self.r * (N - i) *dt)
            FpN[M][0]= 0 
            Fp[:,[i]]= np.dot(inv(M1),FpN)

        #4
        k= np.int32(np.floor(self.S/dS))
        #5
        C = Fc[k][0] + (Fc[k+1][0]-Fc[k][0])/dS * (self.S - k * dS)
        P = Fp[k][0] + (Fp[k+1][0]-Fp[k][0])/dS * (self.S - k * dS)
        return (C, P)