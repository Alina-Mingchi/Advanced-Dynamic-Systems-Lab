#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: felicity
Mingchi Hou
ADSL hw ode
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


############################Problem_1##########################################

#ODE solver
class ODESolver(object):
#Initialize the object with the given initial values
    def __init__(self, p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300):
        self.p_0 = p_0
        self.q_0 = q_0
        self.h = h
        self.n_iter = n_iter

#Explicit Euler method
    def ex_euler(self):
        self.time_ = np.zeros(self.n_iter)
        self.p_ = np.zeros(self.n_iter)
        self.q_ = np.zeros(self.n_iter)
        self.p_[0] = self.p_0
        self.q_[0] = self.q_0
        
        for i in range(self.n_iter-1):
            self.time_[i+1] = self.time_[i] + self.h
            self.p_[i+1] = self.p_[i] + self.h*(-np.sin(self.q_[i]))
            self.q_[i+1] = self.q_[i] + self.h*self.p_[i]
        return self
 
#Implicit Euler method
    def im_euler(self):
        self.time_ = np.zeros(self.n_iter)
        self.p_ = np.zeros(self.n_iter)
        self.q_ = np.zeros(self.n_iter)
        self.p_[0] = self.p_0
        self.q_[0] = self.q_0
        
        for i in range(self.n_iter-1):
            self.time_[i+1] = self.time_[i] + self.h
            self.p_[i+1] = self.p_[i] + self.h*(-np.sin(self.q_[i]))
            self.q_[i+1] = self.q_[i] + self.h*self.p_[i+1]
        return self
    
#Implicit Midpoint method    
    def im_midpoint(self):
        self.time_ = np.zeros(self.n_iter)
        self.p_ = np.zeros(self.n_iter)
        self.q_ = np.zeros(self.n_iter)
        self.p_[0] = self.p_0
        self.q_[0] = self.q_0
        
        for i in range(self.n_iter-1):
            self.time_[i+1] = self.time_[i] + self.h
            self.p_[i+1] = self.p_[i] + self.h*(-np.sin(self.q_[i]))
            self.q_[i+1] = self.q_[i] + 0.5*self.h*(self.p_[i]+self.p_[i+1])
        return self
    
    


#Explicit Euler method
time1=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).ex_euler().time_
q1=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).ex_euler().q_
p1=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).ex_euler().p_
l1 = plt.plot(time1,q1,lw=2,color='blue',label='Explicit Euler')


#Implicit Euler method
time2=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).im_euler().time_
q2=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).im_euler().q_
p2=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).im_euler().p_
l2 = plt.plot(time2,q2,lw=2,color='red',label='Implicit Euler')


#Implicit Midpoint method
time3=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).im_midpoint().time_
q3=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).im_midpoint().q_
p3=ODESolver(p_0 = 1.9, q_0 = 0, h=0.1, n_iter=300).im_midpoint().p_
l3 = plt.plot(time3,q3,lw=2,color='green',label='Implicit Midpoint')


#SciPy build in solver
def f(t, r):
    p = r[0]
    q = r[1]
    return np.array([-np.sin(q), p])

time = np.linspace(0, 30, 1000)
init_r = [1.9, 0]
results = solve_ivp(f, (0, 30), init_r, method='RK45', t_eval=time, rtol=1e-8)
l4 = plt.plot(results.t, results.y[1],lw=2,color='pink',label='Build in solver')



plt.legend()
plt.show()

print("From the above comparison, the Implicit Euler method seems to be closer to the build in solver")
print("The Implicit Midpoint method has a better result compared to the Explicit Euler solver, but both are not stable")



#############################Problem_2#########################################
#Explicit Euler method

plt.figure()
ll1 = plt.plot(q1,p1,lw=2,color='blue',label='Explicit Euler trajectory')
ll2 = plt.plot(q2,p2,lw=2,color='red',label='Implicit Euler trajectory')
ll3 = plt.plot(q3,p3,lw=2,color='green',label='Implicit Midpoint trajectory')

qq = linspace(-10,5,25)
pp = linspace(-5,5,25)
qqq,ppp = np.meshgrid(qq,pp)
EE = 0.5 * np.multiply(ppp,ppp) - np.cos(qqq)
plt.contourf(qqq,ppp,EE)


plt.legend()
plt.show()


print("From the contour and trajectory plot, we can see that the trajectory is surrounding the contour")


############################Problem_3##########################################

#T = 10, build in solver
ttime = np.linspace(0, 10, 1000)
iinit_r = [1.9, 0]
rresults = solve_ivp(f, (0, 10), iinit_r, method='RK45', t_eval=ttime, rtol=1e-8)

plt.figure()
plt.plot(rresults.t, rresults.y[1],lw=2,color='pink',label='Build in solver')
plt.show()








############################Problem_4##########################################

#SciPy build in solver
def f(t, r):
    y = r[0]
    x = r[1]
    mu = 1000
    return np.array([mu*(1-x*x)*y-x, y])

time = np.linspace(0, 10, 1000)
init_r = [1.9, 0]
results = solve_ivp(f, (0, 10), init_r, method='RK45', t_eval=time, rtol=1e-8)
plt.figure()
plt.plot(results.t, results.y[1],lw=2,color='pink',label='Build in solver')
plt.show()


###############################################################################
#References: http://user.astro.columbia.edu/~gbryan/W4260_2010/W4260_lec5.pdf
#https://medium.com/modern-physics/simple-pendulum-odesolver-using-python-dcb30c267eee
#https://stackoverflow.com/questions/56153628/using-scipys-solve-ivp-to-solve-non-linear-pendulum-motion
























