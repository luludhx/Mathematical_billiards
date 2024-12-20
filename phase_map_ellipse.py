# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 14:13:52 2023

@author: dhume
"""

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(14,10))

# Semi major/minor axis
a,b = (2,1)

# Initial conditions: (initial positions, direction wrt horizontal)
thetas_phis = [(np.arctan2(0.5,-np.sqrt(3)/a), np.pi/2),(np.arctan2(0.5,np.sqrt(3)/a),np.pi/2), (np.pi/2, 3*np.pi/2), (0,np.pi), (np.pi,0), (2*np.pi/3, 4*np.pi/3), (2*np.pi/3, 5*np.pi/12), (2*np.pi/3, 7*np.pi/4), (2*np.pi/3, 5*np.pi/3), (2*np.pi/3, 3*np.pi/2), (np.pi, 7*np.pi/12), \
               (np.pi, -7*np.pi/12), (np.pi, 2*np.pi/3), (np.pi, -2*np.pi/3), (np.pi,3*np.pi/4), (np.pi,-3*np.pi/4), (np.pi, 5*np.pi/6), (np.pi, -5*np.pi/6), (np.pi, 11*np.pi/12), (np.pi, -11*np.pi/12)]

# Generating the trajectories for each IC
for theta_phi in thetas_phis:
    theta, phi = theta_phi
    pos = (a * np.cos(theta), b * np.sin(theta))
    
    # Number of trajectories
    steps = 200
    
    Phi = np.zeros((steps,))
    Theta = np.zeros((steps,))
    counter = 0
    while counter < steps:
        Phi[counter] = phi
        Theta[counter] = np.arctan2(pos[1], pos[0])
        
        # Incoming velocity
        v_in = np.array([np.cos(phi), np.sin(phi)])
        
        # Distance to travel for next intersection
        distance = -2 * (pos[0] * v_in[0] / a**2 + pos[1] * v_in[1] / b**2) * (1 / ((v_in[0] / a)**2 + (v_in[1] / b)**2))
        
        # Updating the position
        pos = pos + distance * v_in
        
        # Correcting errors
        e = (pos[0] / a)**2 + (pos[1] / b)**2-1
        if abs(e) < 1e-16:
            x = np.sign(pos[0]) * a * np.sqrt(1 - (pos[1]**2/b**2))
            y = np.sign(pos[1]) * b * np.sqrt(1 - (pos[0]**2 / a**2))
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
        counter +=1
    
    # Plotting the Poincaré section for each IC
    ax.scatter(Phi[1:-1], Phi[2:], marker='.')
plt.show()