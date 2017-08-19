import numpy as np
import matplotlib.pyplot as plt
import math 
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def bs_implicit_exp():
        print()
        print('Black-Scholes PDE, Numerical Estimation, Explicit Finite Differences')
        print()
#Set the current stock price S curr, strike price K, volatility sigma, risk-free rate r, time toexpiration T
        S_curr = 105.99
        K = 105.00
        sigma = 0.3128
        r = 0.0017
        T = 351.0/365.0
        divYield = 0.017
#Set number of asset price steps J
        J = 200
#Set number of time steps M and time increment deltaT
        M = 15000
#Set size of deltaS and deltaT
        S_max = 4*S_curr
        deltaS = S_max/J
        deltaT = T/M
        
        print('Strike Price K = %s' %K)
        print('Volatility sigma = %s' %sigma)
        print('Risk-Free Rate r = %s' %r)
        print('Time to expiration T = %s' %T)
        print('Dividend yield divYeild = %s' %divYield)
        print()
        
        print('Asset price steps J = %s' %J)
        print('Time increment steps M = %s' %M)
        print('Size of Asset price steps = %s' %deltaS)
        print('Size of the time steps = %s' %deltaT)
        
        
#Create (J-1) X (J-1) Matrix A 
#Create column vector of asset prices S
#Create column vector of option prices C_hat
        B=np.zeros((J-1,J-1),dtype=float)
        S=np.zeros((J-1,1),dtype=float)
        C_hat_call=np.zeros((J-1,1),dtype=float)
        C_hat_put=np.zeros((J-1,1),dtype=float)
        
        S[0] = deltaS
         
        j=0 #row counter
        
        for j in range(1,J-1):
            S[j] = S[j-1] + deltaS
        
        print('S = %s' %S)
        print()
        
        
        j=0 #row counter
        for j in range(0,J-1):
            C_hat_call[j] = max(S[j]-K,0)
            C_hat_put[j] = max(K-S[j],0)
            
        print('C_hat_call = %s' %C_hat_call)
        print('C_hat_put = %s' %C_hat_put)
        print()
            
        j=0 #row counter
        k=0 #column counter
        for j in range (0,J-1):
            for k in range (0,J-1):
                if k == j:
                    B[j,j] = 1 + (sigma**2)*((j+1)**2)*deltaT + r*deltaT
                elif k == j-1:
                    B[j,k] = -0.5 * ((sigma**2)*((j+1)**2)*deltaT - (r - divYield)*(j+1)*deltaT)
                elif k == j+1:
                    B[j,k] = -0.5 * ((sigma**2)*((j+1)**2)*deltaT + (r - divYield)*(j+1)*deltaT)
            
        print('B= %s' %B)
        print()
        
        Binverse = np.linalg.matrix_power(B,-1)
        print('Binverse= %s' %Binverse)
        C_0 = 0
        C_J = 0
        C_start_call = C_hat_call
        C_start_put = C_hat_put
        m=1 #time increment counter
        while m <= M:
        #Boundary values calculated separately
        #Must incorporate lowest and highest option prices, which are not in matrix
            C_min_call = 0.5*((sigma**2)*((1)**2)*deltaT - (r - divYield)*(1)*deltaT)*C_0 + Binverse[0,0]*C_hat_call[0] + Binverse[0,1]*C_hat_call[1]
            C_max_call = Binverse[J-2,J-3]*C_hat_call[J-3] + Binverse[J-2,J-2]*C_hat_call[J-2] + 0.5*((sigma**2)*((J-1)**2)*deltaT + (r - divYield)*(J-1)*deltaT)*(S_max - K*math.exp(-r*m*deltaT))
        
            C_min_put = 0.5*((sigma**2)*((1)**2)*deltaT - (r - divYield)*(1)*deltaT)*(K*math.exp(-r*m*deltaT)) + Binverse[0,0]*C_hat_put[0] + Binverse[0,1]*C_hat_put[1]
            C_max_put = Binverse[J-2,J-3]* C_hat_put[J-3] + Binverse[J-2,J-2]*C_hat_put[J-2]  + 0.5*((sigma**2)*((J-1)**2)*deltaT + (r - divYield)*(J-1)*deltaT)*C_J
        #Perform matrix multiplication Ac_M = c_m+1
            C_hat_call = Binverse.dot(C_hat_call)
            C_hat_put = Binverse.dot(C_hat_put)
            
        #Update boundary values
            C_hat_call[0] = C_min_call
            C_hat_call[J-2] = C_max_call
            C_hat_put[0] = C_min_put
            C_hat_put[J-2] = C_max_put
                 
        #Capture option prices during backward walk for graphing
            if m == M/4:
                C1_call = C_hat_call
                C1_put = C_hat_put
            elif m == 2*M/4:
                C2_call = C_hat_call
                C2_put = C_hat_put
            elif m == 3*M/4:
                C3_call = C_hat_call
                C3_put = C_hat_put
            elif m == 4*M/4:
                C4_call = C_hat_call
                C4_put = C_hat_put
                
            m=m+1
            
        print("C_hat_call = %s" %C_hat_call)
        print("C_hat_put = %s" %C_hat_put)
        print()
        print("Asset Price = %s" %S[(J-2)/4])
        print('Call Price = %s' %C_hat_call[(J-2)/4])
        print('Put Price = %s' %C_hat_put[(J-2)/4])
        
        plt.plot(S,C_start_call,'b--',S,C1_call,'r--',S,C2_call,'g--',S,C3_call,'r--',S,C4_call,'b--')
        plt.plot(S,C_start_put,'b--',S,C1_put,'r--',S,C2_put,'g--',S,C3_put,'r--',S,C4_put,'b--')
    

# main program starts here
bs_implicit_exp()