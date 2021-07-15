#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BZ Reaction

Mingchi Hou
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

H = 2 #H_RZ84

k1 = 1e7
k2 = 15 / H**2
k3 = 1.7e4 * H
k4 = 100
k_4 = 1.2e5 / 2
k5 = 1.2e5 / H
kj = 2e-5
q = 0.5

s1 = 800
s2 = 100
s3 = 50
s4 = 25
water = 50

total = s1 + s2 + s3 + s4 + water

A = 0.5 * s1 / total
B = s2 / total
C = 0.025 * s4 / total
H = 0.5 * s1 / total

# Define the three ODEs that we need
def bz(ini,t,A,B,C,H,q):
    X = ini[0]
    Y = ini[1]
    Z = ini[2]
    dX = -k1*H*X*Y + k2*H**2*A*Y - 2*k3*X**2 + k4*H*A*X
    dY = -k1*H*X*Y - k2*H**2*A*Y + q*(kj*B/H)*(Z/(C-Z))
    dZ = -(kj*B/H)*(Z/(C-Z)) + 2*k4*H*A*X
    return dX,dY,dZ


# This period function is presented in Class by Professor Oliver
def period(A,B,C,H,q,Y0):
    t = np.linspace(0,time,num)
    ini = [X0,Y0,Z0]
    sol = odeint(bz,ini,t,args=(A,B,C,H,q))
    s = sol[discard:,0]
    smid = (max(s)+min(s))/2
    a = s>smid
    r = np.logical_and(np.logical_not(a[:-1]),a[1:])
    ii = np.nonzero(r)[0]
    return (t[ii[-1]] - t[ii[0]])/(len(ii)-1)

#Initial values of the ODE
#The initial values are the same as the demo in class
X0 = 0
Y0 = s3/total
Z0 = 0

time = 1000
num = 10000
discard = 1000
pts = 5 #number of points

#Plot the concentration over time
plt.figure()
ini = [X0,Y0,Z0]
t = np.linspace(0,time,num)
concentration = odeint(bz,ini,t,args=(A,B,C,H,q))
plt.plot(t,concentration[:,0],label = 'Concentration of X')
plt.plot(t,concentration[:,2],label = 'Concentration of Z')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Concentration over time')



# %% Varying concentration in solution 1

AA = 0.5 * np.array([800,700,600,500,400])/total
PA = [6.6, 21.1, 76.6, 137.1, 182.8] 

#Simulation
A1 = 0.5*np.linspace(400,800,pts)/total
P1 = [period(a,B,C,a,q,Y0) for a in A1]
# A and H are all from solution 1

plt.figure()
plt.plot(A1,P1,'r-',label = 'Simulation Result')
plt.plot(AA,PA,'ro',label = 'Experiment Result')
plt.legend()
plt.title('Varying Solution 1')
plt.xlabel('Concentration')
plt.ylabel('Period')

# %% Varying concentration in solution 2

BB = np.array([150, 125, 100, 75, 50])/total
PB = [12.9, 8.8, 8.8, 6.6, 31]

B2 = np.linspace(50,150,pts)/total
P2 = [period(A,b,C,H,q,Y0) for b in B2]

plt.figure()
plt.plot(B2,P2,'g-',label = 'Simulation Result')
plt.plot(BB,PB, 'go', label = 'Experiment Result')
plt.legend()
plt.title('Varying Solution 2')
plt.xlabel('Concentration')
plt.ylabel('Period')

# %% Varying concentration in solution 3

CC = np.array([100, 80, 60, 40, 20])/total
PC= [76.4, 38.5, 7.9, 5.8, 19.8]

C3 = np.linspace(20,100,pts)/total
P3 = [period(A,B,C,H,q,y) for y in C3]

plt.figure()
plt.plot(C3,P3,'b-',label = 'Simulation Result')
plt.plot(CC,PC, 'bo', label = 'Experiment Result')
plt.legend()
plt.title('Varying Solution 3')
plt.xlabel('Concentration')
plt.ylabel('Period')


































