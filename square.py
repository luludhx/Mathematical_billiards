# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 14:41:24 2023

@author: dhume
"""
import numpy as np
import matplotlib.pyplot as plt

# Initial position
pos = np.array([0,0.5])
# Direction of trajectory wrt horizontal
phi = np.pi/4

#ang = np.pi/7
# Number of trajectories
steps = 4

def intercepts(pos, v_in):
    '''
    
    Parameters
    ----------
    pos : position on the boundary
    ang : direction wrt horizontal

    Returns
    -------
    A tuple of distances t, such that pos + vt = boundary, for boundaries:
        0. x=0, y=y
        1. x=x, y=0
        2. x=1, y=y
        3. x=x, y=1

    ''' 
    x = pos[0]
    y = pos[1]
    t_0 = -x / v_in[0]
    t_1 = -y / v_in[1]
    t_2 = (1 - x) / v_in[0]
    t_3 = (1 - y) / v_in[1]   
    return t_0, t_1, t_2, t_3

# Normals
n = np.array([(-1,0), (0,-1), (1,0), (0,1)])
# Tangents
t = np.array([(0,1), (1,0), (0,1), (1,0)])

def PosUpdate(distances, pos, v_in, tol = 1e-10):
    '''
    Parameters
    ----------
    distances : array of distances to travel from pos in direction v_in before 
        hitting a boundary
    pos : current position
    v_in : current velocity
    tol : Error margin. The default is 1e-12

    Returns
    -------
    pos : updated position
    index : keeps track of which boundary the current position is on

    '''
    for index, distance in enumerate(distances):
        
        # Getting all the possible positions
        p = pos + distance * v_in
        
        # Checking which position is in the box to 1e-12 accuracy (0<x,y<1)
        if -tol <= p[0] <= 1+tol and -tol<= p[1] <=1+tol and abs(distance) > tol:   
            pos = p
            break
    return pos, index
        
Pos = np.zeros((steps,2))
Phi = np.zeros(steps)
counter = 0
while counter < steps:
    Pos[counter,:] = pos[:]
    Phi[counter] = phi
    
    # Incoming velocity
    v_in = np.array([np.cos(phi), np.sin(phi)])
    
    distances = intercepts(pos, v_in)
    
    pos, index = PosUpdate(distances, pos, v_in)
            
    # Outgoing velocity
    v_out = -np.dot(v_in, n[index]) * n[index] + np.dot(v_in, t[index]) * t[index]
    
    # New direction given by velocity
    phi = np.arctan2(v_out[1], v_out[0])
    
    counter +=1
    
# Plotting the box
fig, ax = plt.subplots(1, 2, figsize=(20,20))
ax[0].plot(np.ones((50,)), np.linspace(0,1,50), color='black')
ax[0].plot(np.zeros((50,)), np.linspace(0,1,50), color='black')
ax[0].plot(np.linspace(0,1,50), np.ones((50,)), color='black')
ax[0].plot(np.linspace(0,1,50), np.zeros((50,)), color='black')

# Plotting trajectories
ax[0].plot(Pos[:,0],Pos[:,1], color='purple')
ax[1].scatter(Phi[1:-1], Phi[2:], marker='.', color='purple')
ax[0].set_aspect('equal', 'box')
ax[1].set_aspect('equal', 'box')
ax[1].set_visible(False)
plt.show()















