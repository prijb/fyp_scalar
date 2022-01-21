#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 00:48:00 2022

@author: prijith
"""

import numpy as np
import numpy.fft as f
import scipy.integrate as integrate
import matplotlib.pyplot as plt
#%%

x0 = 0          #Gaussian mean position
x_sig = 1       #Gaussian standard dev in x
p_sig = 1/x_sig #Gaussian standard dev in p

x = np.linspace(-10,10,100)
t = 0           #Time at which rho and j are calculated

'''
Additional code to try to make it a 2d colour plot
t = np.linspace(-10,10,10)

xx,tt = np.meshgrid(x,t)
'''

#Momentum integration range if resorting to FFTs
px = np.linspace(-10*p_sig,10*p_sig,1000)
py = np.linspace(-10*p_sig,10*p_sig,1000)

'''
Bunch of functions defined as below:

renorm(x,y): Renormalises y(x) over range x
gaussian(x,x0,sig): Defines a normalised gaussian with mean x0 and std sig
integrand_rho: Integrand which has to be double integrated over px and py for rho 
integrand_j: Integrand which has to be double integrated over px and py for j

Note: Currently using psi(x) = sqrt(gaussian(x,x0,x_sig))
                      psi(p) ~ gaussian(p,p0,1/(sqrt(2)*sig))
'''
def renorm(x,y):
    norm = integrate.simpson(y,x)
    y = y/norm
    return y

def gaussian(x,x0,sig):
    y = (1/np.sqrt(2*np.pi*sig*sig))*np.exp(-(x-x0)**2/(2*sig*sig))
    return y

def integrand_rho(px,py,p0,psig,x,t):
    i = (1 + 1/(np.sqrt(1+px**2)*np.sqrt(1+py**2)) + px*py/(np.sqrt(1+px**2)*np.sqrt(1+py**2)))*(np.cos((px-py)*x) - ((np.sqrt(1+px**2)-np.sqrt(1+py**2))*t))*np.exp(-(px**2+py**2)/psig**2)
    
    return i

def integrand_j(px,py,p0,psig,x,t):
    j = (px/np.sqrt(1+px**2) + py/np.sqrt(1+py**2))*(np.cos((px-py)*x - ((np.sqrt(1+px**2)-np.sqrt(1+py**2))*t)))*np.exp(-(px**2+py**2)/psig**2)
    
    return j    

'''
Double integrates the integrands above over px,py ranges [-10,10]
'''
def rho(x,t):
    dens = integrate.nquad(integrand_rho,[[-10,10],[-10,10]],args=(0,p_sig,x,t))
    return dens[0]

def j(x,t):
    dens = integrate.nquad(integrand_j,[[-10,10],[-10,10]],args=(0,p_sig,x,t))
    return dens[0]

'''
Code block for grids (incomplete)
z = np.zeros_like(xx)

for k,t_in in enumerate(t):
    for i,x_in in enumerate(x):
        z[i,k] = rho(x_in,t_in) - np.abs(j(x_in,t_in))
'''

#Code for filling rho and j arrays with x
p_dens = np.zeros_like(x)
j_dens = np.zeros_like(x)

for i,x_in in enumerate(x):
    p_dens[i] = rho(x_in,t)
    j_dens[i] = j(x_in,t)
    
#Rescale rho and j to make sure they're normalised
p_dens = p_dens/integrate.simpson(p_dens,x)
j_dens = j_dens/integrate.simpson(p_dens,x)
#%%
#Plot them against each other
plt.plot(x,p_dens,color='red',label='Probability density') 
plt.plot(x,j_dens,color='blue',label='Probability current density')
plt.xlabel('Position')
plt.ylabel(r'$\rho$/j')
plt.grid()
plt.legend()
plt.show()   
#%%
#Plot for contour 
plt.contourf(xx,tt,z)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.show()

