'''
Mingchi Hou 
code for exercise from lecture 1
'''

import numpy as np
import matplotlib.pyplot as plt


#Initialization of constants
TimeStep = 1        # in years
waterDepth = 4000   # in meters
L = 1350            # solar luminosity in W / m2
alpha = 0.3         # albedo
epsilon = 1         # emissivity
sigma = 5.67e-8     # Stefan-Boltzmann constant in W / m2 K4
s_y = 3.14e7        # seconds / year
HeatCapacity = waterDepth * 4.2e6        # J / K m2 

iteration = 5000

Var = np.zeros(iteration+1) #Variation of heat content
HeatContent = np.zeros(iteration+1)
T = np.zeros(iteration+1)
T0 = np.zeros(iteration+1)
T300 = np.zeros(iteration+1)


###Initial T is 0
T0[0] = 0 
HeatContent[0] = T0[0] * HeatCapacity

for i in range(iteration):
    Var[i] = L * (1 - alpha) / 4 - epsilon * sigma * T0[i] ** 4
    HeatContent[i+1] = HeatContent[i] + Var[i] * TimeStep * s_y
    T0[i+1] = HeatContent[i+1] / HeatCapacity



###Initial T is 300
T300[0] = 300 
HeatContent[0] = T300[0] * HeatCapacity

for i in range(iteration):
    Var[i] = L * (1 - alpha) / 4 - epsilon * sigma * T300[i] ** 4
    HeatContent[i+1] = HeatContent[i] + Var[i] * TimeStep * s_y
    T300[i+1] = HeatContent[i+1] / HeatCapacity



###Initial T is 500
T[0] = 500 
HeatContent[0] = T[0] * HeatCapacity

for i in range(iteration):
    Var[i] = L * (1 - alpha) / 4 - epsilon * sigma * T[i] ** 4
    HeatContent[i+1] = HeatContent[i] + Var[i] * TimeStep * s_y
    T[i+1] = HeatContent[i+1] / HeatCapacity

T500 = T

#Plot
Time = np.linspace(0,iteration,num = iteration+1)  
plt.plot(Time,T0,color='blue',label = "Initial with 0")
plt.plot(Time,T300,color='green',label = "Initial with 300")
plt.plot(Time,T500,color='red',label = "Initial with 500")
plt.xlabel('Time (years)')
plt.ylabel('Planetary temperature (K)')
plt.title('Time dependent temperature')
plt.legend()
plt.show()


































