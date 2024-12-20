# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:22:50 2023

@author: dhume
"""

import numpy as np
import matplotlib.pyplot as plt


# Semi major/minor axis
a,b = (2, 1)
delta = 1e-8

# Initial position
theta = np.pi
pos = (a * np.cos(theta), b * np.sin(theta))
# Direction of trajectory wrt horizontal
phi = delta
# Number of trajectories
steps = 1000

Pos = np.zeros((steps,2))
Phi = np.zeros((steps,))
counter = 0
while counter < steps:
    Pos[counter,:] = pos[:]
    Phi[counter] = phi
    
    # Incoming velocity
    v_in = np.array([np.cos(phi), np.sin(phi)])
   
    # Distance to travel for next intersection
    distance = -2 * (pos[0] * v_in[0] / a**2 + pos[1] * v_in[1] / b**2) * (1 / ((v_in[0] / a)**2 + (v_in[1] / b)**2))
    
    # Updating the position
    pos = pos + distance * v_in
    
    # Correcting errors
    e = (pos[0] / a)**2 + (pos[1] / b)**2 - 1
    
    if abs(e) < 1e-16:
        x = np.sign(pos[0]) * a * np.sqrt(1 - (pos[1] / b)**2)
        y = np.sign(pos[1]) * b * np.sqrt(1 - (pos[0] / a)**2)
        pos = np.array([x, y])
    
    # Normal
    n = np.array([pos[0] / a**2, pos[1] / b**2])
    # Tangent
    t = np.array([pos[1] / b**2, -pos[0] / a**2])
    # Normalising
    n = n / np.linalg.norm(n) 
    t = t / np.linalg.norm(t)
    
    # Outgoing velocity
    v_out = -np.dot(v_in, n) * n + np.dot(v_in, t) * t
    
    phi = np.arctan2(v_out[1], v_out[0])
    
    counter += 1
    
# Plotting the ellipse
fig, ax = plt.subplots(1, 2, figsize=(24,20))
ang = np.linspace(0, 2*np.pi, 100000)
ax[0].plot(a * np.cos(ang), b * np.sin(ang), color='black')
ax[0].set_aspect('equal', 'box')
ax[0].set_xlim(-2.5, 2.5)
ax[0].set_ylim(-1.5, 1.5)

# Distance from the center of the ellipse to each focus
c = np.sqrt(a**2 - b**2)
# Plotting the foci
ax[0].scatter([c,-c], [0,0], color='black', marker='x')

# Plotting the trajectories
ax[0].plot(Pos[:,0], Pos[:,1], color='purple')
# ax[0].set_visible(False)

# Plotting the phase map- previous angle against current
ax[1].scatter(Phi[1:-1], Phi[2:], marker='.', color='purple')
ax[1].set_aspect('equal', 'box')
ax[1].set_visible(False)

plt.show

plt.show()