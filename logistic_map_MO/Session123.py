#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: felicity
Mingchi Hou
ADSL hw logistic map
"""

import numpy as np
import matplotlib.pyplot as plt


############################Problem_1##########################################

def logistic(r, x):
    return r * x * (1 - x)

m = 500 #number of intervals in [0,1] 
n = 10000 #number of chosen r
r = np.linspace(2.0, 4.0, n)

iterations = 1000 #number of iterations
last = 100 #only taking the result from the last 100 iterations

x = 0.5 * np.ones(n)

fig, ax1 = plt.subplots(1, 1, figsize=(8, 8),sharex=True)

#plot the bifurcation diagram
for i in range(iterations):
    x = logistic(r, x)

    if i >= (iterations - last):
        ax1.plot(r, x, ',k', alpha=.05)
ax1.set_xlim(2, 4)
ax1.set_title("Bifurcation diagram")
ax1.set_xlabel('r')
ax1.set_ylabel('x')
plt.savefig('bifurcation.png')

#############################Problem_2#########################################

iterations = 1000 #number of iterations
last = 300 #only taking the result from the last 100 iterations

rr = 3.569945672
x0 = 0.5
x_values = np.zeros(iterations+1) 
x_values[0] = x0



#Compute a matrix of the sequence of x 
for i in range(iterations):
    x1 = rr * x0 * (1-x0)
    x_values[i+1] = x1
    x0 = x1


x_values = x_values[-last:]
            
N_bins = 2**np.arange(1,15,1) #number of bins
epsi = 1/N_bins #size of bin

N_epsi = np.zeros(len(N_bins))  #number of nonzero bins

for j in range(len(N_bins)):
    interval = np.linspace(0.0,1.0,N_bins[j]+1)

    ls = np.histogram(x_values, interval)
    arr = np.matrix.flatten(ls[0])
    N_epsi[j] = np.count_nonzero(arr)


#calculate linear regression
coeff = np.polyfit(np.log(1/epsi),np.log(N_epsi),1)    

#Plot the capacity dimension
fig, ptt = plt.subplots(1, 1, figsize=(8,8))
ptt1 = ptt.plot(np.log(1/epsi),np.log(N_epsi),'o',label="Actual data") 
ptt2 = ptt.plot(np.log(1/epsi),np.polyval(coeff,np.log(1/epsi)),'r',label="Linear regression")
   
ptt.set_title("Capacity Dimension")
ptt.set_xlabel('-log(epsi)')
ptt.set_ylabel('log(N_epsi)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
plt.show()
fig.savefig('dimension.png')

print("The capacity dimention is the slope from the plot, we can also extract from the linear regression")
print(str(coeff[0]) + " is the capacity dimension")


############################Problem_3##########################################

rrr = np.linspace(3.5, 3.7, 1000)
x0 = 0.5

d = np.zeros(len(rrr)) 

for k in range(len(rrr)):

    r_temp = rrr[k]
    #Compute a matrix of the sequence of x 
    x_val = np.zeros(iterations+1) 
    x_val[0] = x0
    for i in range(iterations):
        x1 = r_temp * x0 * (1-x0)
        x_val[i+1] = x1
        x0 = x1
            
        N_bins = 2**np.arange(1,15,1) #number of bins
        epsi = 1/N_bins #size of bin

    N_epsi = np.zeros(len(N_bins))  #number of nonzero bins
    x_val = x_val[-last:]
    for j in range(len(N_bins)):
        interval = np.linspace(0.0,1.0,N_bins[j]+1)

        ls = np.histogram(x_val, interval)
        arr = np.matrix.flatten(ls[0])
        N_epsi[j] = np.count_nonzero(arr)


#calculate linear regression
    coeff = np.polyfit(np.log(1/epsi),np.log(N_epsi),1)    
    d[k] = coeff[0]
    
    
###Plot the capacity dimension
fig, cdr = plt.subplots(1, 1, figsize=(8,8))
cdr.plot(rrr,d,'o')
cdr.set_title("Capacity Dimension around critical value")
cdr.set_xlabel('r')
cdr.set_ylabel('capacity dimension')
plt.savefig('dimension_r.png')

print("I choose to use maximal 16384 bins, and thus smallest bin size is 6.10351562e-05")
print("There are 1000 iterations, but I use the last 300 ones, i.e.700 preiterations")


print("We need to choose a small enough bin size so as to approximate the limit better, while the smaller the bin size gets, the longer it takes to compute")
print("With too few preiterations, the capacity dimension tends to be larger than the actual")
print("With not enough iterations, the capacity dimension tends to be smaller than the actual")
###############################################################################
#Reference: https://stackoverflow.com/questions/2946519/counting-number-of-values-between-interval
#Reference: https://ipython-books.github.io/121-plotting-the-bifurcation-diagram-of-a-chaotic-dynamical-system/


















