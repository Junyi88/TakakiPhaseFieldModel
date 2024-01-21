#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 16:48:53 2024

@author: junyi
"""

import numpy as np
import matplotlib.pyplot as plt

#%%

dtheta = np.linspace(0, 2.0*np.pi/0.001, num=10001)
beta = 1.0e5
omega = 0.0016
mu = 600.0

q = 1.0 - np.exp(-beta * omega * dtheta) + (mu / omega) * np.exp(-beta*omega*dtheta)
q1 = np.exp(-beta * omega * dtheta)
q2 = (mu / omega) * np.exp(-beta*omega*dtheta)

plt.figure(1, clear=True)
plt.plot(dtheta, q, 'b-')
plt.plot(dtheta, q1, 'r-')
plt.plot(dtheta, q2, 'k-')


plt.figure(2, clear=True)
plt.plot(dtheta, 3.0e-5*1.0e6/q, 'b-')
# plt.plot(dtheta, 1.0/q1, 'r-')
# plt.plot(dtheta, 1.0/q2, 'k-')
