#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diffusion 1

Mingchi Hou
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

def analytical(x):
    return - x ** 2 / 2 + 1 / 2

def L2error(array,xarr):
    temp = array[:-1]
    term1 = (array[1:] + temp) / 2
    xtemp = xarr[:-1]
    xtempi = xarr[1:]
    term2 = analytical((xtemp+xtempi)/2)
    E0 = np.sqrt((2/len(xtemp))*np.sum(np.abs(term1-term2) ** 2))
    return E0

def Lerror(array,xarr):
    temp = array[:-1]
    term1 = (array[1:] + temp) / 2
    xtemp = xarr[:-1]
    xtempi = xarr[1:]
    term2 = analytical((xtemp+xtempi)/2)
    E1 = max(np.abs(term1-term2))
    return E1
    
    

# Define constants
g = 0
a = 2
f = 2

j = np.array([2,3,4,5,6,7,8,9,10])
N = 2 ** j
h = 2 / N

E0arr = np.zeros(len(j))
E1arr = np.zeros(len(j))

# part c
# Compute analytical solution and numerical solution 
for num in j:
    print(num)

    farr = np.zeros(N[num-2]+1) #\vec{f}
    # Set boundary conditions
    farr[0] = g
    farr[-1] = g
    
    for i in range(N[num-2]-1):
        farr[i+1] = f
    
    
    M = np.zeros((N[num-2]+1,N[num-2]+1))
    M[0,0] = 1
    M[-1,-1] = 1
    for i in range(N[num-2]-1):
        M[i+1,i] = 1
        M[i+1,i+1] = -2
        M[i+1,i+2] = 1
    M =  a / h[num-2] ** 2 *M
    
    
    Minv = inv(M)
    uarr = - Minv.dot(farr)
    
    x_an = np.arange(-1,1+0.01,0.01) # Analytical
    x_nu = np.arange(-1,1+h[num-2],h[num-2]) #Numerical
    
    plt.figure()
    plt.plot(x_an,analytical(x_an),label = 'Analytical Result')
    plt.plot(x_nu,uarr,label = 'Numerical Result')
    plt.legend()
    plt.title('Plot j = ' + str(num))
    plt.xlabel('x')
    plt.ylabel('u(x)')


# part d
    #Calculate errors
    E0 = L2error(uarr,x_nu)
    # print(E0)
    E0arr[num-2] = E0
    
    E1 = Lerror(uarr,x_nu)
    # print(E1)
    E1arr[num-2] = E1

plt.figure()
plt.yscale('log')
plt.plot(j,E0arr,label = 'log E0')
plt.plot(j,E1arr,label = 'log E1')
plt.legend()
plt.title('Log error plot')
plt.xlabel('j')
plt.ylabel('Error')

logE0arr = np.log(E0arr)
logE1arr = np.log(E1arr)

slopeE0 = (logE0arr[-1]-logE0arr[0])/8
# print(E0arr)
# print(logE0arr)
print('The slope of E0 is' + str(slopeE0))
slopeE1 = (logE1arr[-1]-logE1arr[0])/8
# print(E1arr)
# print(logE1arr)
print('The slope of E1 is' + str(slopeE1))


# %% part e

def ana(x):
    ana_value = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0:
            ana_value[i] =  - x[i] ** 2  - x[i]/2 + 1/2
        else:
            ana_value[i] =  - x[i]/2 + 1/2
    return ana_value

def L2error(array,xarr):
    temp = array[:-1]
    term1 = (array[1:] + temp) / 2
    xtemp = xarr[:-1]
    xtempi = xarr[1:]
    term2 = ana((xtemp+xtempi)/2)
    E0 = np.sqrt((2/len(xtemp))*np.sum(np.abs(term1-term2) ** 2))
    return E0

def Lerror(array,xarr):
    temp = array[:-1]
    term1 = (array[1:] + temp) / 2
    xtemp = xarr[:-1]
    xtempi = xarr[1:]
    term2 = ana((xtemp+xtempi)/2)
    E1 = max(np.abs(term1-term2))
    return E1

# Define constants
g = 0
a = 1

j = np.array([2,3,4,5,6,7,8,9,10])
N = 2 ** j
h = 2 / N

E0arr = np.zeros(len(j))
E1arr = np.zeros(len(j))

for num in j:
    print(num)

    farr = np.zeros(N[num-2]+1) #\vec{f}
    # Set boundary conditions
    farr[0] = g
    farr[-1] = g
    
    for i in range(N[num-2]-1):
        if i < N[num-2]/2 -1:
            farr[i+1] = 2
        else:
            farr[i+1] = 0
    
    print(farr)
    
    M = np.zeros((N[num-2]+1,N[num-2]+1))
    M[0,0] = 1
    M[-1,-1] = 1
    for i in range(N[num-2]-1):
        M[i+1,i] = 1
        M[i+1,i+1] = -2
        M[i+1,i+2] = 1
    M =  a / h[num-2] ** 2 *M
    
    
    Minv = inv(M)
    uarr = - Minv.dot(farr)
    
    x_an = np.arange(-1,1+0.01,0.01) # Analytical
    x_nu = np.arange(-1,1+h[num-2],h[num-2]) #Numerical
    
    plt.figure()
    plt.plot(x_an,ana(x_an),label = 'Analytical Result')
    plt.plot(x_nu,uarr,label = 'Numerical Result')
    plt.legend()
    plt.title('Plot j = ' + str(num))
    plt.xlabel('x')
    plt.ylabel('u(x)')


    #Calculate errors
    E0 = L2error(uarr,x_nu)
    # print(E0)
    E0arr[num-2] = E0
    
    E1 = Lerror(uarr,x_nu)
    # print(E1)
    E1arr[num-2] = E1

plt.figure()
plt.yscale('log')
plt.plot(j,E0arr,label = 'log E0')
plt.plot(j,E1arr,label = 'log E1')
plt.legend()
plt.title('Log error plot')
plt.xlabel('j')
plt.ylabel('Error')

logE0arr = np.log(E0arr)
logE1arr = np.log(E1arr)

slopeE0 = (logE0arr[-1]-logE0arr[0])/8
# print(E0arr)
# print(logE0arr)
print('The slope of E0 is' + str(slopeE0))
slopeE1 = (logE1arr[-1]-logE1arr[0])/8
# print(E1arr)
# print(logE1arr)
print('The slope of E1 is' + str(slopeE1))








