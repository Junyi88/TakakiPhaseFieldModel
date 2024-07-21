# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 20:44:26 2024

@author: Lenovo
"""

import numpy as np
import matplotlib.pyplot as plt

#%%
NY = 11
NX = 1001
x = np.linspace(-0.01, 0.09, 1001)
x = np.linspace(-0.01, 0.09, 1001)
# x = np.linspace(-50.0, 50.0, 1001)
dx = x[1] - x[0]

k = 0.5e3
c = 0.06

theta = 0.2 + 0.25 * (np.tanh(k * (x - c)) + 1.0)
# y = np.tanh(x)

#%%
plt.figure(1, clear=True)
plt.plot(x, theta, 'b-')

#%%
# a = 0.5
mu = 0.06
s = 20.0e-4

f = (1.0 / (s * np.sqrt(2.0 * np.pi)) ) * np.exp(-0.5 * (x - mu) * (x - mu) / (s * s))
f = f / np.max(f)
g = 1.0 - 0.1 * f
plt.plot(x, g/2.0, 'r-')

plt.figure(2, clear=True)
plt.plot(x, f, 'b-')
plt.plot(x, g, 'r-')

#%%
Theta = np.zeros([NY, NX])
Phi = np.zeros([NY, NX])

for i in range(NY):
    Theta[i,:] = theta
    Phi[i,:] = g

np.savetxt("./ignoreDir/theta0.csv", Theta, delimiter=',')
np.savetxt("./ignoreDir/phi0.csv", Phi, delimiter=',')
