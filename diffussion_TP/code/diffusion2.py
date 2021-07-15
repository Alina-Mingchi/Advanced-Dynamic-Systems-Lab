#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diffusion 2

Mingchi Hou
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv


# %% part b
# exp and approximations

value = np.arange(-3,1+0.1,0.1)
EE = 1 + value
IE = 1 / (1-value)
CN = (1+value/2) / (1-value/2)
 
plt.figure()
plt.plot(value,np.exp(value),'black',label = 'exp')
plt.plot(value,EE,'r',label='Explicit Euler')
plt.plot(value,IE,'g',label='Implicit Euler')
plt.plot(value,CN,'b',label='Crank-Nicholson')
plt.legend()
plt.xlim([-3,1])
plt.ylim([-2,3])
plt.title('Plot of actual value and approximations')
plt.xlabel('x')
plt.ylabel('exp(x)')


# %% part c
# Plot with stability regions

# Explicit Euler
plt.figure()
plt.plot(value,np.exp(value),'black',label = 'exp')
plt.plot(value,EE,'r',label='Explicit Euler')
plt.plot(value[np.abs(EE)<1],EE[np.abs(EE)<1],'r*', label='stability region')
plt.legend()
plt.xlim([-3,1])
plt.ylim([-2,3])
plt.title('Explicit Euler and stability region')
plt.xlabel('x')
plt.ylabel('exp(x)')

# Implicit Euler
plt.figure()
plt.plot(value,np.exp(value),'black',label = 'exp')
plt.plot(value,IE,'g',label='Implicit Euler')
plt.plot(value[np.abs(IE)<1],IE[np.abs(IE)<1],'g*', label='stability region')
plt.legend()
plt.xlim([-3,1])
plt.ylim([-2,3])
plt.title('Implicit Euler and stability region')
plt.xlabel('x')
plt.ylabel('exp(x)')

# Crank-Nicholson
plt.figure()
plt.plot(value,np.exp(value),'black',label = 'exp')
plt.plot(value,CN,'b',label='Crank-Nicholson')
plt.plot(value[np.abs(CN)<1],CN[np.abs(CN)<1],'b*', label='stability region')
plt.legend()
plt.xlim([-3,1])
plt.ylim([-2,3])
plt.title('Crank-Nicholson and stability region')
plt.xlabel('x')
plt.ylabel('exp(x)')


# %%part d

def ini(x):
    return (1-x)*(1-x)*(1+x)*(1+x)

def EE(y):
    return 1+y

def IE(y):
    return inv(1-y)

def CN(y):
    Q = 1+y/2
    R = 1-y/2
    Rinv = inv(R)
    return Rinv.dot(Q)

# Define constants
a = 1
T = 0.5
j = 3
N = 2 ** j
h = 2 / N


l = np.array([1,2,3,4,5,6,7,8])
tau = 1 / (2 ** l)
step = 0.5/tau
t = np.linspace(0,T,N+1)

M = np.zeros((N+1,N+1))
M[0,0] = 1
M[-1,-1] = 1
for i in range(N-1):
    M[i+1,i] = 1
    M[i+1,i+1] = -2
    M[i+1,i+2] = 1
M =  a / h ** 2 *M


uini =  ini(t)

for index in range(len(l)):                                        
    # EE
    temp1 = np.identity(N+1)
    for i in range(int(step[index])):
        ee = EE(tau[index]*M)
        temp1 = ee.dot(temp1)
    eeuarr = temp1.dot(uini)    
    
    # IE
    temp2 = np.identity(N+1)
    for i in range(int(step[index])):
        ie = IE(tau[index]*M)
        temp2 = ie.dot(temp2)
    ieuarr = temp2.dot(uini)  
    
    # CN
    temp3 = np.identity(N+1)
    for i in range(int(step[index])):
        cn = CN(tau[index]*M)
        temp3 = cn.dot(temp3)
    cnuarr = temp3.dot(uini) 
    
    # actual
    temp4 = np.identity(N+1)
    for i in range(int(step[index])):
        real = np.exp(tau[index]*M)
        temp4 = real.dot(temp4)
    realuarr = temp4.dot(uini) 
    
    plt.figure()
    plt.plot(t,eeuarr,'r',label='EE')
    plt.plot(t,ieuarr,'g',label='IE')
    plt.plot(t,cnuarr,'b',label='CN')
    plt.legend()
    # plt.xlim([-3,1])
    # plt.ylim([-2,3])
    plt.title('Plot of tau = '+str(tau[index]))
    plt.xlabel('t')
    plt.ylabel('u')

    print(eeuarr)
    print(ieuarr)
    print(cnuarr)
    print(realuarr)
    print('-------------------------------')


# %% part e

# EE #tau = 0.25
temp1 = np.identity(N+1)
for i in range(int(step[1])):
    ee = EE(tau[1]*M)
    temp1 = ee.dot(temp1)
eeu1 = temp1.dot(uini)    

# CN #tau = 0.25
temp3 = np.identity(N+1)
for i in range(int(step[1])):
    cn = CN(tau[1]*M)
    temp3 = cn.dot(temp3)
cnu1 = temp3.dot(uini)    

plt.figure()
plt.plot(t,eeu1,'r',label='EE')
plt.plot(t,cnu1,'b',label='CN')
plt.legend()
# plt.xlim([-3,1])
# plt.ylim([-2,3])
plt.title('Plot of tau = 0.25')
plt.xlabel('t')
plt.ylabel('u')


# EE #tau = 0.015625
temp1 = np.identity(N+1)
for i in range(int(step[5])):
    ee = EE(tau[5]*M)
    temp1 = ee.dot(temp1)
eeu2 = temp1.dot(uini) 

# CN #tau = 0.015625
temp3 = np.identity(N+1)
for i in range(int(step[5])):
    cn = CN(tau[5]*M)
    temp3 = cn.dot(temp3)
cnu2 = temp3.dot(uini) 

plt.figure()
plt.plot(t,eeu2,'r',label='EE')
plt.plot(t,cnu2,'b',label='CN')
plt.legend()
# plt.xlim([-3,1])
# plt.ylim([-2,3])
plt.title('Plot of tau = 0.015625')
plt.xlabel('t')
plt.ylabel('u')















