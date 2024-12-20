# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:56:50 2023

@author: dhume
"""

import numpy as np
import matplotlib.pyplot as plt

# Half length of rectangle
L = 2
# Radius of semi circle
r = 1

# Initial position
pos = np.array([3,0])
# Direction of trajectory wrt horizontal
phi = -np.pi/3
# Number of trajectories
steps = 200000

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
    
    # Left edge circle
    t21 = -x/v_in[0]
    
    #Bottom line
    t_3 = (-r - pos[1]) / v_in[1] 
    
    return t01, t02, t21, t_1, t_3

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
    for distance in distances:
        # Updating the position
        p = pos + distance * v_in
        if L==0:
            if tol < p[0] <= 1+tol and -1-tol< p[1]<1+tol and abs(distance)>1e-6:
                pos = p
                n = np.array([pos[0], pos[1]]) 
                t = np.array([-pos[1], pos[0]])
                break
            elif abs(p[0])<tol and -1-tol< p[1]<1+tol and abs(distance)>1e-6:
                pos = p
                n = np.array([-1, 0]) 
                t = np.array([0, 1]) 
                break
            else:continue
        else:    
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
                
            # Check for left line
            elif abs(p[0]) < tol and -r-tol <= p[1] <=r+tol:
                pos = p
                pos = np.array([0, p[1]])
                n = np.array([1, 0])
                t = np.array([0, -1])
                break
                
            # Check for top line
            elif -tol <= p[0] <= L+tol and abs(p[1] - r) < tol:
                pos = np.array([p[0], r])
                n = np.array([0, -1])
                t = np.array([1, 0])
                break
            
            # Check for bottom line
            elif -tol <= p[0] <= L+tol and abs(p[1] + r) < tol:
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
    distances = [distance for distance in dists if np.abs(distance) > 1e-11]
    
    pos, n, t = PosUpdate(distances, pos, v_in, L, r)
    
    # Outgoing velocity
    v_out = -np.dot(v_in, n) * n + np.dot(v_in, t) * t
    
    phi = np.arctan2(v_out[1], v_out[0])
    counter += 1
    
# Plotting the stadium
fig, ax = plt.subplots(1, 2, figsize=(25,20))
ax[0].set_aspect('equal', 'box')
x = np.array([L,0,0,L])
y = np.array([-r,-r,r,r])
ax[0].plot(x,y,color = "black")
thetas = np.linspace(-np.pi/2,np.pi/2,100)
x = L + r*np.cos(thetas)
y = r*np.sin(thetas)
ax[0].plot(x,y,color = "black")
    
# Plotting the trajectories
ax[0].plot(Pos[:,0], Pos[:,1], color='purple')
ax[0].set_visible(False)

# Plotting the phase space
ax[1].scatter(Theta, Phi, marker='p', s = 1, color = 'purple', alpha = 0.9)
ax[1].set_aspect('0.5', 'box')
#ax[1].set_visible(False)

    
