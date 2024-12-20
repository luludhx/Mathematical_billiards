# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 22:51:48 2023

@author: dhume
"""

import numpy as np
import matplotlib.pyplot as plt

# Initial position 
theta = np.pi
pos = (np.cos(theta), np.sin(theta))
# Direction of trajectory wrt horizontal
phi = np.pi/3.25724267
# Number of trajectories
steps = 100

Pos = np.zeros((steps,2))
Phi = np.zeros((steps,))
Theta = np.zeros((steps,))
counter = 0
while counter < steps:
    Pos[counter,:] = pos[:]
    Phi[counter] = phi
    Theta[counter] = np.arctan2(pos[1], pos[0])
    
    # Incoming velocity
    v_in = np.array([np.cos(phi), np.sin(phi)])
    
    # Distance to travel for next intersection
    distance = -2 * (pos[0] * v_in[0] + pos[1] * v_in[1])
    
    # Updating the position
    pos = pos + distance * v_in
    
    # Normal
    n = np.array([pos[0], pos[1]]) 
    # Tangent
    t = np.array([-pos[1], pos[0]])
    
    # Outgoing velocity
    v_out = -np.dot(v_in, n) * n + np.dot(v_in, t) * t
    
    # Updating direction wrt horizontal
    phi = np.arctan2(v_out[1], v_out[0]) 
    counter += 1

# Plotting the circle
fig, ax = plt.subplots(figsize=(11, 11))
theta = np.linspace(0, 2 * np.pi, 50)
ax.plot(np.cos(theta), np.sin(theta), color='black')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
plt.gca().set_aspect('equal', adjustable='box')

# Plotting the trajectories
plt.plot(Pos[:,0], Pos[:,1], color='purple')

plt.show()