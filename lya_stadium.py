# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:56:50 2023

@author: dhume
"""

import numpy as np
import matplotlib.pyplot as plt

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
    t_1 = (r - pos[1]) / v_in[1]
    
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
    t_3 = (-r - pos[1]) / v_in[1] 
    
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
        distance = -2 * (pos[0] * v_in[0] + pos[1] * v_in[1])
        pos = pos + distance * v_in
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
                    x = -np.sqrt(r**2 - pos[1]**2) - L
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

def Main(pos, phi, steps, L, r):
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
        
    return Pos, Phi, Theta

delta = 1e-12
steps = 200
L = 40
r = 1
PosStart = np.array([-r-L/2,0])
PhiStart = 2*np.pi/3.123456789

Pos, VelAng, PosAng = Main(PosStart, PhiStart, steps, L, r)

PosD, VelAngD, PosAngD = Main(PosStart+delta, PhiStart+delta, steps, L, r)

# The phi Lyapunov
s = np.log(np.sqrt((VelAng - VelAngD)**2 + (PosAng - PosAngD)**2))
a = np.polyfit(np.arange(1,steps+1), s, 1)
print(f'The Lyapunov for phi is {a[0]}')

# The Euclidean distance Lyapunov
s2 = np.log(np.sqrt((Pos[:,0] - PosD[:,0])**2 + (Pos[:,1] - PosD[:,1])**2))
a2 = np.polyfit(np.arange(1,steps+1), s2, 1)
print(f'The Lyapunov for the Euclidean distance is {a2[0]}')

# Plotting the phase space Lyapunov
fig, ax = plt.subplots(1, 5, figsize=(60,65))
ax[0].plot(s)
ax[0].set_aspect('equal', 'box')
ax[0].set_visible(False)

# Plotting the Euclidean Lyapunov
ax[1].plot(s2)
ax[1].set_aspect('equal', 'box')
ax[1].set_visible(False)

# Plotting the angle difference
ax[2].plot(VelAng)
ax[2].plot(VelAngD)
ax[2].set_visible(False)
ax[2].set_aspect('equal', 'box')

# Plotting the distance difference
ax[3].plot(Pos[:,0])
ax[3].plot(PosD[:,0])
ax[3].set_visible(False)
ax[3].set_aspect('equal', 'box')

# A list of Lyapunov for different L found by manually choosing the trajectories (L, lambda)
Lya = np.array([(0, 0.053), (1, 1.076), (3, 0.833), (5,0.741), (10,0.712), (20,0.230), (40, 0.152), (80,0.073)])

#Plotting the Lyapunov against the length L
ax[4].plot(Lya[:,0], Lya[:,1], marker='o', color='purple')
ax[4].set_aspect('25', 'box')
