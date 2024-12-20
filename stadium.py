# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 19:00:36 2023

@author: dhume
"""

import numpy as np
import matplotlib.pyplot as plt

# Half length of rectangle
L = 1
# Radius of semi circles
r = 1

# Initial position
pos = np.array([0,-1])
# Direction of trajectory wrt horizontal
phi = 5*np.pi/6
# Number of trajectories
steps = 100000

def intercepts(pos, v_in, r, L):
    '''
    Parameters
    ----------
    pos : position on the boundary
    v_in : incoming velocity 
    r : Tadius of semi circles
    L : Half length of rectangle

    Returns
    -------
    A tuple of distances t, such that pos + vt = boundary, for boundaries:
        0. (x-L)^2 +y^2 -r^2 =0
        1. x=x, y=r
        2. (x+L)^2 +y^2 -r^2 =0
        3. x=x, y=-r

    ''' 
    x = pos[0]
    y = pos[1]
    
    # Right semi circle
    A0 = 1
    B0 = (2 * (x - L) * v_in[0] + 2 * y * v_in[1])
    C0 = (x - L)**2 - r**2 + y**2
    sr0= B0**2 - 4 * A0 * C0
    if sr0 < 0:
        sr0 = complex(sr0)
    t01 = (-B0 + np.sqrt(sr0)) / 2 * A0
    t02 = (-B0 - np.sqrt(sr0)) / 2 * A0
    
    # Top line
    if v_in[1] != 0:
        t_1 = (r - pos[1]) / v_in[1]
    else: t_1 = 0
    
    # Left semi circle
    A2 = 1
    B2 = 2 * (x + L) * v_in[0] + 2 * y * v_in[1]
    C2 = (x + L)**2 - r**2 + y**2  
    sr2 = B2**2 - 4 * A2 * C2
    if sr2 < 0:
        sr2 = complex(sr2)
    t21 = (-B2 + np.sqrt(sr2)) / 2 * A2
    t22 = (-B2 - np.sqrt(sr2)) / 2 * A2
    
    #Bottom line
    if v_in[1] != 0:
        t_3 = (-r - pos[1]) / v_in[1]
    else: t_3 = 0
    
    return t01, t02, t21, t22, t_1, t_3

def is_complex(x):
    '''
    Parameters
    ----------
    x : scalar

    Returns
    -------
    True if x is a complex number

    '''
    return isinstance(x, complex)

def PosUpdate(distances, pos, v_in, L, r, tol = 1e-10):
    '''
    Parameters
    ----------
    distances : array of distances to travel from pos in direction v_in before 
        hitting a boundary
    pos : current position
    v_in : current velocity
    tol : Error margin. The default is 1e-10

    Returns
    -------
    pos : updated position
    n : normal of boundary on which the updated position is
    t : tangent of boundary on which the updated position is

    '''
    if L == 0:
        d = -2 * (pos[0] * v_in[0] + pos[1] * v_in[1])
        pos = pos + d * v_in
        n = np.array([pos[0], pos[1]]) 
        t = np.array([-pos[1], pos[0]])
        
    else: 
        for distance in distances:
            # Updating the position
            p = pos + distance * v_in
            
            # Check for right semi circle
            if L-tol <= p[0] <= L + r +tol and -r-tol <= p[1] <=r+tol:
                if abs((p[0] - L)**2 + p[1]**2 - r) < tol:
                    pos = p
                    x = np.sqrt(r**2 - pos[1]**2) + L
                    y = np.sign(pos[1]) * np.sqrt(r**2 - (pos[0] - L)**2)
                    pos = np.array([x, y])
                    n = np.array([pos[0] - L, pos[1]])
                    t = np.array([-pos[1], pos[0] - L])
                    break
                
            # Check for left semi circle
            elif -L-r-tol <= p[0] <= -L+tol and -r-tol <= p[1] <=r+tol:
                if abs((p[0] + L)**2 + p[1]**2 - r) < tol:
                    pos = p
                    x = -np.sqrt(abs(r**2 - pos[1]**2)) - L
                    y = np.sign(pos[1]) * np.sqrt(r**2 - (pos[0] + L)**2)
                    pos = np.array([x, y])
                    n = np.array([pos[0] + L, pos[1]])
                    t = np.array([-pos[1], pos[0] + L])
                    break
                
            # Check for top line
            elif -L-tol <= p[0] <= L+tol and abs(p[1] - r) < tol:
                pos = np.array([p[0], r])
                n = np.array([0, 1])
                t = np.array([-1, 0])
                break
            
            # Check for bottom line
            elif -L-tol <= p[0] <= L+tol and abs(p[1] + r) < tol:
                pos = np.array([p[0], -r])
                n = np.array([0, -1])
                t = np.array([1, 0])
                break
            else: continue  
    
    return pos, n, t

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
    
    d = intercepts(pos, v_in, r, L)
    # Getting rid of complex distances
    dists = [distance for distance in d if is_complex(distance) is False]
    # Getting rid of the 0 distance
    distances = [distance for distance in dists if np.abs(distance) > 1e-6]
    
    pos, n, t = PosUpdate(distances, pos, v_in, L, r)
    
    # Outgoing velocity
    v_out = -np.dot(v_in, n) * n + np.dot(v_in, t) * t
    
    phi = np.arctan2(v_out[1], v_out[0])
    counter += 1


# Plotting the stadium
fig, ax = plt.subplots(1, 2, figsize=(25,20))
ax[0].set_aspect('equal', 'box')
# ax[0].set_xlim(-3.3,3.3)
# ax[0].set_ylim(-1.7,1.7)

# Top semi circle
theta_t = np.linspace(-np.pi / 2, np.pi / 2, 1000)
x_t = r * np.cos(theta_t) + L
y_t = r * np.sin(theta_t)
ax[0].plot(x_t, y_t, color='black')

# Bottom semi circle
theta_b = np.linspace(-np.pi / 2, np.pi / 2, 1000)
x_b = r * np.cos(theta_b) + L
y_b = r * np.sin(theta_b)
ax[0].plot(-x_b, y_b, color='black')

# Straight lines
x_s = np.linspace(-L , L , 1000)
y_st = np.full_like(x_s, r)
y_bt = np.full_like(x_s, -r)
ax[0].plot(x_s, y_st, color='black')
ax[0].plot(x_s, y_bt, color='black')

# Plotting the trajectories
ax[0].plot(Pos[:,0], Pos[:,1], color='purple')
ax[0].set_visible(False)

# Plotting the phase space
ax[1].scatter(Theta, Phi, marker='p', s = 1, color = 'purple', alpha = 0.9)
ax[1].set_aspect('equal', 'box')
#ax[1].set_visible(False)

plt.show()







