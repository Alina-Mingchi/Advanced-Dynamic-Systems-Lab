#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mingchi Hou 
Graded exercise 
Modelling daisyworld
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


'''Part 1'''
## Set parameters
# aw,ab area of daisies
# ax area of empty ground
# bw,bb growth rate of daisies
# t time
# Tw,Tb local temperature
# Te global temperature
# Ap albedo of the planet

Ab = 0.25           #Albedo
Ax = 0.5            #Albedo
Aw = 0.75           #Albedo
Topt = 22.5         #optimal temperature
r = 0.3             #death rate of daisies
sigma = 5.6704e-8   #Stefan-Boltzmann constant
S = 917             #solar constant
q = 20              #constant


## Set initial conditions for ODEs, time array 
L = 1.2             #Luminosity
dt = 1
tmax = 100
y0 = [0.1,0.1]
t = np.arange(0,tmax,dt)
dydt = np.zeros((2,tmax))

def func(y,t):
    Ap = Ab*y[1] + Ax*(1-y[0]-y[1]) + Aw*y[0]
    Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15
    Tw = q*(Ap-Aw) + Te
    Tb = q*(Ap-Ab) + Te
    bw = max(0,1 - ((Topt-Tw)/17.5) ** 2)
    bb = max(0,1 - ((Topt-Tb)/17.5) ** 2)
    dy1dt = y[0]*((1-y[0]-y[1])*bw - r) #aw
    dy2dt = y[1]*((1-y[0]-y[1])*bb - r) #ab
    dydt = [dy1dt,dy2dt]
    return dydt

# Solve the ODEs
y = odeint(func,y0,t)        

empty = np.ones((1,tmax))
empty = empty - y[:,0] - y[:,1]
empty = empty.T  #Areas covered by the empty ground

Ap = Ab*y[:,1] + Ax*(1-y[:,0]-y[:,1]) + Aw*y[:,0]
Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15

# Plot
fig,ax = plt.subplots(2)
fig.suptitle('Lunimosity = 1.2')
ax[0].plot(t,y[:,0],'r-',label='white')
ax[0].plot(t,y[:,1],'b-',label='black')
ax[0].plot(t,empty,'g-',label='empty')
ax[0].set_ylim([0,1])
ax[0].set_xlim([0,100])
ax[0].set_ylabel('Area covered')
ax[0].legend(prop={"size":7})
ax[1].plot(t,Te,label='Global Temperature')
ax[1].set_xlim([0,100])
ax[1].set_ylabel('Temperature')
ax[1].set_xlabel('Time')
ax[1].legend()
plt.show()



'''Part 2'''
## Set initial conditions for ODEs, time array 
del L
dt = 1
tmax = 100
y0 = [0.1,0.1]
t = np.arange(0,tmax,dt)
dydt = np.zeros((2,tmax))

def func(y,t):
    Ap = Ab*y[1] + Ax*(1-y[0]-y[1]) + Aw*y[0]
    Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15
    Tw = q*(Ap-Aw) + Te
    Tb = q*(Ap-Ab) + Te
    bw = max(0,1 - ((Topt-Tw)/17.5) ** 2)
    bb = max(0,1 - ((Topt-Tb)/17.5) ** 2)
    dy1dt = y[0]*((1-y[0]-y[1])*bw - r) #aw
    dy2dt = y[1]*((1-y[0]-y[1])*bb - r) #ab
    dydt = [dy1dt,dy2dt]
    return dydt

#Initialize arraies to be plotted
plaw = np.zeros((0))
plab = np.zeros((0))
plax = np.zeros((0))
plTe = np.zeros((0))

#Loop over increasing Luminosity
for L in np.arange(0.6, 1.41, 0.02):

    y = odeint(func,y0,t)  
    Ap = Ab*y[:,1] + Ax*(1-y[:,0]-y[:,1]) + Aw*y[:,0]
    Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15
    
    plaw = np.append(plaw,y[-1,0])
    plab = np.append(plab,y[-1,1])
    plax = np.append(plax,1 - y[-1,0] - y[-1,1])    
    plTe = np.append(plTe,Te[-1])

Larr = np.arange(0.6, 1.41, 0.02)
Tee = (S*Larr*(1-Ax)/sigma) ** 0.25 - 273.15    
#Tee: global temperature no daisy

# Plot
fig2,ax2 = plt.subplots(2)
fig2.suptitle('Plot of part 2')
ax2[0].plot(Larr,plaw,'r-',label='white')
ax2[0].plot(Larr,plab,color='black',label='black')
ax2[0].plot(Larr,plax,'g-',label='empty')
ax2[0].set_ylim([0,1])
ax2[0].set_xlim([0.6,1.4])
ax2[0].set_ylabel('Area covered')
ax2[0].legend(prop={"size":7})
ax2[1].plot(Larr,plTe,label='Global Temperature')
ax2[1].plot(Larr,Tee,color='black',linestyle='dashed',label='Global Temperature (no daisies)')
ax2[1].set_ylim([-10,50])
ax2[1].set_xlim([0.6,1.4])
ax2[1].set_ylabel('Temperature')
ax2[1].set_xlabel('Luminosity')
ax2[1].legend()
plt.show()
#


'''Part 3'''

## Set parameters
# a_i area of species
# ax area of empty ground
# b_i growth rate of species
# t time
# T_i local temperature
# Te global temperature
# Ap albedo of the planet

A1 = 0.1            #Albedo
A2 = 0.2
A3 = 0.3
A4 = 0.4
Ax = 0.5            #Albedo for empty ground
A6 = 0.6
A7 = 0.7
A8 = 0.8
A9 = 0.9         
Topt = 22.5         #optimal temperature
r = 0.3             #death rate of daisies
sigma = 5.6704e-8   #Stefan-Boltzmann constant
S = 917             #solar constant
q = 20              #constant


## Set initial conditions for ODEs, time array 
L = 1.2             #Luminosity
dt = 1
tmax = 1000
y0 = [0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
t = np.arange(0,tmax,dt)
dydt = np.zeros((8,tmax))

def func(y,t):
    ax = 1-y[0]-y[1]-y[2]-y[3]-y[4]-y[5]-y[6]-y[7]
    Ap = A1*y[0] + A2*y[1] + A3*y[2] + A4*y[3]+ A6*y[4] + A7*y[5] + A8*y[6] + A9*y[7]+ Ax*ax 
    Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15
    
    T1 = q*(Ap-A1) + Te
    T2 = q*(Ap-A2) + Te
    T3 = q*(Ap-A3) + Te
    T4 = q*(Ap-A4) + Te
    T6 = q*(Ap-A6) + Te
    T7 = q*(Ap-A7) + Te
    T8 = q*(Ap-A8) + Te
    T9 = q*(Ap-A9) + Te
    
    b1 = max(0,1 - ((Topt-T1)/17.5) ** 2)
    b2 = max(0,1 - ((Topt-T2)/17.5) ** 2)
    b3 = max(0,1 - ((Topt-T3)/17.5) ** 2)
    b4 = max(0,1 - ((Topt-T4)/17.5) ** 2)    
    b6 = max(0,1 - ((Topt-T6)/17.5) ** 2)
    b7 = max(0,1 - ((Topt-T7)/17.5) ** 2)
    b8 = max(0,1 - ((Topt-T8)/17.5) ** 2)
    b9 = max(0,1 - ((Topt-T9)/17.5) ** 2)

    
    dy1dt = y[0]*(ax*b1 - r) #a1
    dy2dt = y[1]*(ax*b2 - r)
    dy3dt = y[2]*(ax*b3 - r) 
    dy4dt = y[3]*(ax*b4 - r)    
    dy6dt = y[4]*(ax*b6 - r) 
    dy7dt = y[5]*(ax*b7 - r)    
    dy8dt = y[6]*(ax*b8 - r) 
    dy9dt = y[7]*(ax*b9 - r)    
    
    dydt = [dy1dt,dy2dt,dy3dt,dy4dt,dy6dt,dy7dt,dy8dt,dy9dt]
    return dydt

# Solve the ODEs
y = odeint(func,y0,t)        

empty = np.ones((1,tmax))
empty = empty - y[:,0] - y[:,1] - y[:,2] - y[:,3] - y[:,4] - y[:,5] - y[:,6] - y[:,7]
empty = empty.T  #Areas covered by the empty ground

ax = 1-y[:,0]-y[:,1]-y[:,2]-y[:,3]-y[:,4]-y[:,5]-y[:,6]-y[:,7]
Ap = A1*y[:,0] + A2*y[:,1] + A3*y[:,2] + A4*y[:,3]+ A6*y[:,4] + A7*y[:,5] + A8*y[:,6] + A9*y[:,7]+ Ax*ax 

Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15

# Plot
fig3,ax3 = plt.subplots(2)
fig3.suptitle('Lunimosity = 1.2')
ax3[0].plot(t,y[:,0],label='Species1')
ax3[0].plot(t,y[:,1],label='Species2')
ax3[0].plot(t,y[:,2],label='Species3')
ax3[0].plot(t,y[:,3],label='Species4')
ax3[0].plot(t,y[:,4],label='Species5')
ax3[0].plot(t,y[:,5],label='Species6')
ax3[0].plot(t,y[:,6],label='Species7')
ax3[0].plot(t,y[:,7],label='Species8')
ax3[0].plot(t,empty,color='black',label='empty')
ax3[0].set_ylim([0,1])
ax3[0].set_xlim([0,1000])
ax3[0].set_ylabel('Area covered')
ax3[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
ax3[1].plot(t,Te,label='Global Temperature')
ax3[1].set_xlim([0,1000])
ax3[1].set_ylabel('Temperature')
ax3[1].set_xlabel('Time')
ax3[1].legend()
plt.show()

#Plot the Areas coved by differnt species over time
#The same color as compared to the previous plot 
#shows the same species
plt.plot(y)
plt.title('Plot of different Species')
plt.ylabel('Area covered')
plt.xlabel('Time')
plt.show()





'''Optional'''
del L #Delete the definition of L
del y 
dt = 1
tmax = 800
y0 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
t = np.arange(0,tmax,dt)
dydt = np.zeros((2,tmax))


def func(y,t):
    ax = 1-y[0]-y[1]-y[2]-y[3]-y[4]-y[5]-y[6]-y[7]
    Ap = A1*y[0] + A2*y[1] + A3*y[2] + A4*y[3]+ A6*y[4] + A7*y[5] + A8*y[6] + A9*y[7]+ Ax*ax 
    Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15
    
    T1 = q*(Ap-A1) + Te
    T2 = q*(Ap-A2) + Te
    T3 = q*(Ap-A3) + Te
    T4 = q*(Ap-A4) + Te
    T6 = q*(Ap-A6) + Te
    T7 = q*(Ap-A7) + Te
    T8 = q*(Ap-A8) + Te
    T9 = q*(Ap-A9) + Te
    
    b1 = max(0,1 - ((Topt-T1)/17.5) ** 2)
    b2 = max(0,1 - ((Topt-T2)/17.5) ** 2)
    b3 = max(0,1 - ((Topt-T3)/17.5) ** 2)
    b4 = max(0,1 - ((Topt-T4)/17.5) ** 2)    
    b6 = max(0,1 - ((Topt-T6)/17.5) ** 2)
    b7 = max(0,1 - ((Topt-T7)/17.5) ** 2)
    b8 = max(0,1 - ((Topt-T8)/17.5) ** 2)
    b9 = max(0,1 - ((Topt-T9)/17.5) ** 2)

    
    dy1dt = y[0]*(ax*b1 - r) #a1
    dy2dt = y[1]*(ax*b2 - r)
    dy3dt = y[2]*(ax*b3 - r) 
    dy4dt = y[3]*(ax*b4 - r)    
    dy6dt = y[4]*(ax*b6 - r) 
    dy7dt = y[5]*(ax*b7 - r)    
    dy8dt = y[6]*(ax*b8 - r) 
    dy9dt = y[7]*(ax*b9 - r)    
    
    dydt = [dy1dt,dy2dt,dy3dt,dy4dt,dy6dt,dy7dt,dy8dt,dy9dt]
    return dydt

#Initialize arraies to be plotted
pla1 = np.zeros((0))
pla2 = np.zeros((0))
pla3 = np.zeros((0))
pla4 = np.zeros((0))
pla6 = np.zeros((0))
pla7 = np.zeros((0))
pla8 = np.zeros((0))
pla9 = np.zeros((0))
plax = np.zeros((0))
plTe = np.zeros((0))

#Loop over increasing Luminosity

'uncomment this following two lines to implement increasing luminosity'
#Larr = np.arange(0.6, 1.41, 0.02)
#for L in np.arange(0.6, 1.41, 0.02):
'uncomment this following two lines to implemtne decreasing luminosity'
Larr = np.arange(1.41, 0.6, -0.02)
for L in np.arange(1.41, 0.6, -0.02):    
    y = odeint(func,y0,t)  
    
    
    ax = 1-y[:,0]-y[:,1]-y[:,2]-y[:,3]-y[:,4]-y[:,5]-y[:,6]-y[:,7]
    Ap = A1*y[:,0] + A2*y[:,1] + A3*y[:,2] + A4*y[:,3]+ A6*y[:,4] + A7*y[:,5] + A8*y[:,6] + A9*y[:,7]+ Ax*ax 
    Te = (S*L*(1-Ap)/sigma) ** 0.25 - 273.15
    
    pla1 = np.append(plaw,y[-1,0])
    pla2 = np.append(plab,y[-1,1])
    pla3 = np.append(plaw,y[-1,2])
    pla4 = np.append(plab,y[-1,3])
    pla6 = np.append(plaw,y[-1,4])
    pla7 = np.append(plab,y[-1,5])
    pla8 = np.append(plaw,y[-1,6])
    pla9 = np.append(plab,y[-1,7])    
    
    plax = np.append(plax,ax[-1])    
    plTe = np.append(plTe,Te[-1])


Tee = (S*Larr*(1-Ax)/sigma) ** 0.25 - 273.15    
#Tee: global temperature with no species

# Plot
fig4,ax4 = plt.subplots(2)
#fig4.suptitle('Increasing Luminosity')
ax4[0].plot(Larr,pla1[0:-1],label='Species 1')
ax4[0].plot(Larr,pla2[0:-1],label='Species 2')
ax4[0].plot(Larr,pla3[0:-1],label='Species 3')
ax4[0].plot(Larr,pla4[0:-1],label='Species 4')
ax4[0].plot(Larr,pla6[0:-1],label='Species 5')
ax4[0].plot(Larr,pla7[0:-1],label='Species 6')
ax4[0].plot(Larr,pla8[0:-1],label='Species 7')
ax4[0].plot(Larr,pla9[0:-1],label='Species 8')
ax4[0].plot(Larr,plax,color='black',label='empty')
ax4[0].set_ylim([0,1])
ax4[0].set_xlim([0.6,1.4])
ax4[0].set_ylabel('Area covered')
ax4[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.55),
          ncol=3, fancybox=True, shadow=True)
ax4[1].plot(Larr,plTe,label='Global Temperature')
ax4[1].plot(Larr,Tee,color='black',linestyle='dashed',label='Global Temperature (no daisies)')
#ax4[1].set_ylim([-10,50])
ax4[1].set_xlim([0.6,1.4])
ax4[1].set_ylabel('Temperature')
ax4[1].set_xlabel('Luminosity')
ax4[1].legend()
plt.show()







########Reference############
#https://apmonitor.com/pdc/index.php/Main/SolveDifferentialEquations
