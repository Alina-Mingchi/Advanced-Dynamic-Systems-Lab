#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diffusion 3

Mingchi Hou
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

N = 5
N = 25
h = 1/N

M = np.zeros(((N+1)**2,(N+1)**2))
M[0,0] = -4
M[0,1] = 1
M[0,N] = 1

for i in range(N-1):
    M[i+1,i] = 1
    M[i+1,i+1] = -4
    M[i+1,i+2] = 1
    M[i+1,N+i+1] = 1

for ii in range(N**2+1):
    print(ii+N)
    M[ii+N,ii] = 1
    M[ii+N,ii+N-1] = 1
    M[ii+N,ii+N] = -4 
    M[ii+N,ii+N+1] = 1
    M[ii+N,ii+N+N] = 1

for iii in range(N-1):
    M[-(iii+2),-(iii+1)] = 1
    M[-(iii+2),-(iii+2)] = -4
    M[-(iii+2),-(iii+3)] = 1
    M[-(iii+2),-(N+iii+2)] = 1

M[-1,-1] = -4
M[-1,-2] = 1
M[-1,-(N+1)] = 1

M =  1 / h ** 2 *M
Minv = inv(M)

farr = np.zeros((N+1)**2)

for y in range(N+1):
    for x in range(N+1):
        farr[x+y*(N+1)] = 2*(1-x/5)*x/5 + 2*(1-y/5)*y/5

uarr = - Minv.dot(farr)


plt.figure()
plt.plot(uarr)
plt.legend()
plt.title('Plot of u')
plt.xlabel('index')
plt.ylabel('u')




# %% Bunus part e
N = 5
e = 0.5
e = 0.8
e = 0.2
h = 1/N

M = np.zeros(((N+1)**2,(N+1)**2))
M[0,0] = -2-2*e
M[0,1] = 1
M[0,N] = 1

for i in range(N-1):
    M[i+1,i] = 1
    M[i+1,i+1] = -2-2*e
    M[i+1,i+2] = 1
    M[i+1,N+i+1] = 1

for ii in range(N**2+1):
    print(ii+N)
    M[ii+N,ii] = 1
    M[ii+N,ii+N-1] = 1
    M[ii+N,ii+N] = -2-2*e
    M[ii+N,ii+N+1] = 1
    M[ii+N,ii+N+N] = 1

for iii in range(N-1):
    M[-(iii+2),-(iii+1)] = 1
    M[-(iii+2),-(iii+2)] = -2-2*e
    M[-(iii+2),-(iii+3)] = 1
    M[-(iii+2),-(N+iii+2)] = 1

M[-1,-1] = -2-2*e
M[-1,-2] = 1
M[-1,-(N+1)] = 1

M =  1 / h ** 2 *M
Minv = inv(M)

farr = np.zeros((N+1)**2)

for y in range(N+1):
    for x in range(N+1):
        farr[x+y*(N+1)] = 2*(1-x/5)*x/5 + 2*(1-y/5)*y/5

uarr = - Minv.dot(farr)


plt.figure()
plt.plot(uarr)
plt.legend()
plt.title('Plot of u')
plt.xlabel('index')
plt.ylabel('u')




