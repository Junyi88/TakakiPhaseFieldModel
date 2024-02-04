#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 10:22:01 2024

@author: junyi
"""

import numpy as np
import matplotlib.pyplot as plt
import os

#%%
step_start = 10
step_end = 100000

res_dir = "/home/junyi/scratch/diffusionbonding/run3/"
run_name = "Test3" 

theta_start_path = os.path.join(res_dir, "{}_Theta_{}.csv".format(run_name, step_start))
theta_end_path = os.path.join(res_dir, "{}_Theta_{}.csv".format(run_name, step_end))

phi_start_path = os.path.join(res_dir, "{}_Phi_{}.csv".format(run_name, step_start))
phi_end_path = os.path.join(res_dir, "{}_Phi_{}.csv".format(run_name, step_end))

#%%
phi_start = np.loadtxt(phi_start_path, delimiter=',')
theta_start = np.loadtxt(theta_start_path, delimiter=',')

phi_end = np.loadtxt(phi_end_path, delimiter=',')
theta_end = np.loadtxt(theta_end_path, delimiter=',')

#%%
plt.figure(101, clear=True)
phi_0 = phi_start[0,:]
plt.plot(phi_0, 'r-')
phi_1 = phi_end[0,:]
plt.plot(phi_1, 'b--')

plt.figure(102, clear=True)
theta_0 = theta_start[0,:]
plt.plot(theta_0, 'r-')
theta_1 = theta_end[0,:]
plt.plot(theta_1, 'b--')