# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 14:26:19 2023

@author: schec
"""

import numpy as np
import matplotlib.pyplot as plt

def Radius(theta,a=1,b=1):
    return  1/np.sqrt(np.cos(theta)**2 / a**2 +  np.sin(theta)**2 / b**2)

def PlotBunny(delta = 0):
    #head
    x1 = np.linspace(-1.5125613407560135,1.5125613407560135,500)
    y1 = -(x1**2 /2) +3.2
    plt.plot(x1,y1, color='black')
    
    x2 = np.linspace(-1.0682444417944663, 1.0682444417944663,500)
    y2 = (x2**2/2.5)-1.5
    plt.plot(x2,y2, color='black')
    
    #left ear
    y3 = np.linspace(-3.454163456597992, 0.3, 500)
    x3 = (-y3**2/15)-0.1*y3-1.1
    plt.plot(x3,y3, color='black')
    
    y4 = np.linspace(-3.454163456597992, 2.056079095225186, 500)
    x4 = (y4**2/15)+0.1*y4-2
    plt.plot(x4,y4, color='black')
    
    #right ear
    y5 = np.linspace(-3.454163456597992, 0.3, 50)
    x5 = (y5**2/15)+0.1*y5+1.1
    plt.plot(x5,y5, color='black')
    
    y6 = np.linspace(-3.454163456597992, 2.056079095225186, 500)
    x6 = (-y6**2/15)-0.1*y6+2
    plt.plot(x6,y6, color='black')
    
    # Left Eye
    theta = np.linspace(0,2*np.pi,1000)
    x = Radius(theta,0.18,0.23)*np.cos(theta) - 0.6
    y = Radius(theta,0.18,0.23)*np.sin(theta) + 1.1
    plt.plot(x,y, color = "black")
    
    # Right Eye
    theta = np.linspace(0,2*np.pi,1000)
    x = Radius(theta,0.18,0.23)*np.cos(theta) + 0.6
    y = Radius(theta,0.18,0.23)*np.sin(theta) + 1.1
    plt.plot(x,y, color = "black")
    
    
    plt.fill(0.6+ 0.18*np.cos(theta),1.08+ 0.23*np.sin(theta), color='black')
    plt.fill(-0.6+ 0.18*np.cos(theta),1.08+ 0.23*np.sin(theta), color='black')
    
    #Nose Translation
    
    plt.plot([0,-0.2,0.2,0],[-0.2+delta,0.1+delta,0.1+delta,-0.2+delta], color='black')
    plt.xlim(-4,4)
    plt.ylim(-4.5,4.5)


# %% Outline
fig, ax = plt.subplots(figsize=(10,9))

#head
x1 = np.linspace(-1.5125613407560135,1.5125613407560135,500)
y1 = -(x1**2 /2) +3.2
ax.plot(x1,y1, color='black')

x2 = np.linspace(-1.0682444417944663, 1.0682444417944663,500)
y2 = (x2**2/2.5)-1.5
ax.plot(x2,y2, color='black')

#left ear
y3 = np.linspace(-3.454163456597992, 0.3, 500)
x3 = (-y3**2/15)-0.1*y3-1.1
ax.plot(x3,y3, color='black')

y4 = np.linspace(-3.454163456597992, 2.056079095225186, 500)
x4 = (y4**2/15)+0.1*y4-2
ax.plot(x4,y4, color='black')

#right ear
y5 = np.linspace(-3.454163456597992, 0.3, 50)
x5 = (y5**2/15)+0.1*y5+1.1
ax.plot(x5,y5, color='black')

y6 = np.linspace(-3.454163456597992, 2.056079095225186, 500)
x6 = (-y6**2/15)-0.1*y6+2
ax.plot(x6,y6, color='black')

# Left Eye
theta = np.linspace(0,2*np.pi,1000)
x = Radius(theta,0.18,0.23)*np.cos(theta) - 0.6
y = Radius(theta,0.18,0.23)*np.sin(theta) + 1.1
ax.plot(x,y, color = "black")

# Right Eye
theta = np.linspace(0,2*np.pi,1000)
x = Radius(theta,0.18,0.23)*np.cos(theta) + 0.6
y = Radius(theta,0.18,0.23)*np.sin(theta) + 1.1
ax.plot(x,y, color = "black")

'''
# Circle Mouth?
theta = np.linspace(0,2*np.pi,500)
x = Radius(theta,0.3,0.2)*np.cos(theta)
y = Radius(theta,0.3,0.2)*np.sin(theta) - 0.8
plt.plot(x,y, color = "black")
'''

ax.fill(0.6+ 0.18*np.cos(theta),1.08+ 0.23*np.sin(theta), color='black')
ax.fill(-0.6+ 0.18*np.cos(theta),1.08+ 0.23*np.sin(theta), color='black')

ax.plot([0,-0.2,0.2,0],[-0.2,0.1,0.1,-0.2], color='black')
ax.set_xlim(-4,4)
ax.set_ylim(-4.5,4.5)

# %% Plot
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

def intercepts(pos, v_in, delta = 0):
    '''
   
    Parameters
    ----------
    pos : position on the boundary
    ang : direction wrt horizontal

    Returns
    -------
    A tuple of distances t, such that pos + vt = boundary, for boundaries:
        1. y1 = -(x1**2 /2) +3.2
            [-1.736,1.736]
        2. y2 = (x2**2/2.5)-1.5
            [-1.101, 1.101]
        3. x3 = (-y3**2/15)-0.1*y3-1.1
            [-3.719, 0.3]
        4. x4 = (y4**2/15)+0.1*y4-2
            [-3.719, 1.992]
        5. x5 = (y5**2/15)+0.1*y5+1.1
            [-3.719, 0.3]
        6. x6 = (-y6**2/15)-0.1*y6+2
            [-3.719, 1.992]
    '''
    

    
    x=pos[0]
    y=pos[1]
    a=v_in[0]
    b=v_in[1]
   
    A1 = a**2/2
    B1 = x*a+b
    C1 = y+(x**2/2)-3.2
    sr1= B1**2 - 4 * A1 * C1
    if sr1 < 0:
        sr1 = complex(sr1)
    t11 = (-B1 + np.sqrt(sr1)) / (2 * A1)
    t12 = (-B1 - np.sqrt(sr1)) / (2 * A1)

    A2 = -a**2*0.4
    B2 = -x*a*0.8+b
    C2 = y-0.4*x**2+1.5
    sr2= B2**2 - 4 * A2 * C2
    if sr2 < 0:
        sr2 = complex(sr2)
    t21 = (-B2 + np.sqrt(sr2)) / (2 * A2)
    t22 = (-B2 - np.sqrt(sr2)) / (2 * A2)
       
    
    A3 = -b**2/15
    B3 = -2*y*(b/15)-0.1*b-a
    C3 = -(y**2/15)-x-0.1*y-1.1
    sr3= B3**2 - 4 * A3 * C3
    if sr3 < 0:
        sr3 = complex(sr3)
    t31 = (-B3 + np.sqrt(sr3)) / (2 * A3)
    t32 = (-B3 - np.sqrt(sr3)) / (2 * A3)
     

    A4 = b**2/15
    B4 = 2*y*(b/15)+0.1*b-a
    C4 = y**2/15-x+0.1*y-2
    sr4= B4**2 - 4 * A4 * C4
    if sr4 < 0:
        t41 = -10
        t42 = -10
        #sr4 = complex(sr4)
    else:
        t41 = (-B4 + np.sqrt(sr4)) / (2 * A4)
        t42 = (-B4 - np.sqrt(sr4)) / (2 * A4)        
    
    A5 = b**2/15
    B5 = 2*y*(b/15)+0.1*b-a
    C5 = (y**2/15)-x+0.1*y+1.1
    sr5 = B5**2 - 4 * A5 * C5
    if sr5 < 0:
        t51 = -10
        t52 = -10
    else:
        t51 = (-B5 + np.sqrt(sr5)) / (2 * A5)
        t52 = (-B5 - np.sqrt(sr5)) / (2 * A5)  
        
    A6 = -b**2/15
    B6 = -2*y*(b/15)-0.1*b-a
    C6 = -(y**2/15)-x-0.1*y+2
    sr6= B6**2 - 4 * A6 * C6
    if sr6 < 0:
        sr6 = complex(sr6)
    t61 = (-B6 + np.sqrt(sr6)) / (2 * A6)
    t62 = (-B6 - np.sqrt(sr6)) / (2 * A6)
    
    # Left Eye
    k = 1
    A7 = a**2/(0.18)**2 + b**2/(0.23)**2
    B7 = (2*a*x + k*1.2*a)/0.18**2 + (2*b*y - 2.2*b)/0.23**2
    C7 = (0.36 + x**2 + k*1.2*x)/0.18**2 + (1.21 + y**2 - 2.2*y)/0.23**2 -1
    sr7 = B7**2 - 4 * A7 * C7
    if sr7 < 0:
        t71 = -10
        t72 = -10
    else:
        t71 = (-B7 + np.sqrt(sr7)) / (2 * A7)
        t72 = (-B7 - np.sqrt(sr7)) / (2 * A7)
        
    # Right Eye
    k = -1
    A8 = a**2/(0.18)**2 + b**2/(0.23)**2
    B8 = (2*a*x + k*1.2*a)/0.18**2 + (2*b*y - 2.2*b)/0.23**2
    C8 = (0.36 + x**2 + k*1.2*x)/0.18**2 + (1.21 + y**2 - 2.2*y)/0.23**2 -1
    sr8 = B8**2 - 4 * A8 * C8
    if sr8 < 0:
        t81 = -10
        t82 = -10
    else:
        t81 = (-B8 + np.sqrt(sr8)) / (2 * A8)
        t82 = (-B8 - np.sqrt(sr8)) / (2 * A8)
    
    
    # Nose: (1)
    t1 = (0.1 - y - delta) / b
    
    # Nose: (2)
    t2 = (-15*x + 10*(y- delta) + 2)/(15*a - 10*b)
    
    # Nose: (3)
    t3 = -(15*x + 10*(y- delta) + 2)/(15*a + 10*b)
   
    return t11,t12,t21,t22,t31,t32,t41,t42,t51,t52,t61,t62, t71, t72, t81, t82, t1, t2, t3

def Curves(x,y, delta =0):
    C1 = y + x**2 /2 - 3.2
    C2 = y - x**2 / 2.5 + 1.5
    C3 = x + y**2 /15 + y/10 + 1.1
    C4 = x - y**2 /15 -y/10 + 2
    C5 = x - y**2 /15 - y/10 - 1.1
    C6 = x + y**2 /15 + y/10 - 2
    C7 = ((x + 0.6)/0.18)**2 + ((y - 1.1)/0.23)**2-1
    C8 = ((x - 0.6)/0.18)**2 + ((y - 1.1)/0.23)**2-1
    # Nose(1)
    C9 = (y + delta - 0.1)
    # Nose(2)
    C10 = y - delta + 0.2 + 1.5*x
    # Nose(3)
    C11 = y - delta + 0.2 - 1.5*x
    
    
    return C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11

def DrawCurrent(pos,v_in, l = 1):
    s = np.linspace(0,l,100)
    plt.scatter(pos[0],pos[1])
    plt.plot(pos[0] + s*v_in[0],pos[1] + s*v_in[1])
    
    

def PosUpdate(distances, pos, v_in, tol = 1e-7, delta = 0):
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
    val = []
    for distance in distances:
        # Getting all the possible positions
        p = pos + distance * v_in
        x=p[0]
        y=p[1]
        if counter == 15:
            print("hi")
            #s = np.linspace(0,0.2,100)
            #plt.plot(pos[0] + s*v_in[0],pos[1] + s*v_in[1])
            #plt.scatter(pos[0],pos[1])
        # Check 1
        if -1.5125613407560135-tol < x < 1.5125613407560135+tol and 2.05607909-tol < y < 3.2+tol and abs(Curves(x,y)[0]) < tol:
            n = np.array([-x,-1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,1))
        # Check 2
        if -1.0682144417-tol < x < 1.0682144417+tol and -1.5-tol < y < -1.04354152+tol and abs(Curves(x,y)[1]) < tol:
            n = np.array([4*x/5,-1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,2))
        # Check 3
        if -1.55-tol < x < -17/16+tol and -3.454163456-tol < y < 0.3+tol and abs(Curves(x,y)[2]) < tol:
            n = np.array([-1,-(2*y/15)-0.1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,3))
        # Check 4 negative
        if (-2.1 < x < -1.55+tol) and -3.454163456-tol < y < 0 and abs(Curves(x,y)[3]) < tol:
            n = np.array([-1,(2*y/15)+0.1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,4))
            
        # Check 4 positive
        if (-2.1 < x < -1.512613407+tol) and 0 < y < 2.05607909 and abs(Curves(x,y)[3]) < tol:
            n = np.array([-1,(2*y/15)+0.1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,4))
            
        
        # Check 5
        if 1.55+tol > x > -tol + 17/16 and -3.454163456-tol < y < 0.3+tol and abs(Curves(x,y)[4]) < tol:
            n = np.array([-1,(2*y/15)+0.1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,5))
        
        # Check 6 negative
        if (2.1 > x > 1.55-tol) and (-3.454163456-tol < y < 0) and abs(Curves(x,y)[5]) < tol:
            n = np.array([-1,-(2*y/15)-0.1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,6))
            
        # Check 6 positive
        if (2.1 > x > -1.512613407-tol) and 0 < y < 2.05607909 and abs(Curves(x,y)[5]) < tol:
            n = np.array([-1,-(2*y/15)-0.1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,6))
            

        # Check 7
        if (-0.4 > x > -0.8+tol) and 0.8 < y < 1.35 and abs(Curves(x,y)[6]) < tol:
            n = np.array([2*x/0.18**2 + 1.2/0.18**2,2*y/0.23**2 - 2.2/0.23**2])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,7))
        
        # Check 8
        if (0.4 < x < 0.8+tol) and 0.8 < y < 1.35 and abs(Curves(x,y)[7]) < tol:
            n = np.array([2*x/0.18**2 -1.2/0.18**2,2*y/0.23**2 - 2.2/0.23**2])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,8))
        
        # Check 9
        if (-0.2-tol < x < 0.2+tol) and 0.1-tol < (y+delta) < 0.1 + tol and abs(Curves(x,y,delta=delta)[8]) < tol:
            n = np.array([0,-1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,9))

        # Check 10
        if (-0.2-tol < x < tol) and -0.2-tol < (y+delta) < 0.1 + tol and abs(Curves(x,y,delta=delta)[9]) < tol:
            n = np.array([-1.5,-1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,10))
        
        # Check 11
        if (-tol < x < 0.2+tol) and -0.2-tol < (y+delta) < 0.1 + tol and abs(Curves(x,y,delta=delta)[10]) < tol:
            n = np.array([1.5,-1])
            t = np.array((-n[1],n[0]))
            n = n / np.linalg.norm(n)
            t = t / np.linalg.norm(t)
            val.append((p,n,t,distance,11))
        
    
    if len(val)==0:
        print('it fucked up at', counter)
    elif len(val)==1:
        return val[0][0:-2]
    else:
        diffs = []
        for v in val:
            p, _, _, d,_ = v
            diffs.append(d)
        
        val = np.array(val)
        diffs = np.array(diffs)
        val = val[diffs > 0]
        diffs = diffs[diffs > 0]
        if len(diffs) < 1:
            print('it fucked up at', counter)
        index = np.argmin(diffs)
        
        return  val[index][0:-2]
    

pos1 = (0,-1.5)
pos2 = (0,-1.5)

steps = 100
phi1 = np.pi/2.05
dphi = 1e-6
phi2 = phi1 + dphi

#phi = -np.pi/4

T = np.linspace(0,int(steps / 10),steps)

#T = T[0:4]
#steps = 4

Sniff = False
if Sniff:
    delta = 0.02*np.sin(8*T)**2
else:
    delta = np.zeros(steps)
    

Pos1 = np.zeros((steps,2))
Phi1 = np.zeros(steps)
Pos2 = np.zeros((steps,2))
Phi2 = np.zeros(steps)


counter = 0
Animate = True
while counter < steps:
    Pos1[counter,:] = pos1[:]
    Phi1[counter] = phi1
    
    Pos2[counter,:] = pos2[:]
    Phi2[counter] = phi2
    
   
    # Incoming velocity
    v_in1 = np.array([np.cos(phi1), np.sin(phi1)])
    v_in2 = np.array([np.cos(phi2), np.sin(phi2)])
   
    d1 = intercepts(pos1, v_in1, delta = delta[counter])
    d2 = intercepts(pos2, v_in2, delta = delta[counter])
    # Getting rid of complex distances
    dists1 = [distance for distance in d1 if is_complex(distance) is False]
    dists2 = [distance for distance in d2 if is_complex(distance) is False]
    
    # Getting rid of the 0 distance
    distances1 = [distance for distance in dists1 if np.abs(distance) > 1e-7]
    distances2 = [distance for distance in dists2 if np.abs(distance) > 1e-7]
    #distances = dists
    pos1, n1, t1 = PosUpdate(distances1, pos1, v_in1, delta = delta[counter])
    pos2, n2, t2 = PosUpdate(distances2, pos2, v_in2, delta = delta[counter])
   
    # Outgoing velocity
    v_out1 = -np.dot(v_in1, n1) * n1 + np.dot(v_in1, t1) * t1
    v_out2 = -np.dot(v_in2, n2) * n2 + np.dot(v_in2, t2) * t2
   
    phi1 = np.arctan2(v_out1[1], v_out1[0])
    phi2 = np.arctan2(v_out2[1], v_out2[0])
    
    
    if Animate:
        PlotBunny(delta = delta[counter])
        if counter > 1:
            plt.plot(Pos1[0:counter,0],Pos1[0:counter,1],alpha = 0.2, color = "purple")
            plt.plot(Pos1[counter-1:counter+1,0],Pos1[counter-1:counter+1:,1], color='purple')
            
            plt.plot(Pos2[0:counter,0],Pos2[0:counter,1],alpha = 0.2, color = "indigo")
            plt.plot(Pos2[counter-1:counter+1,0],Pos2[counter-1:counter+1:,1], color='indigo')
            
        #plt.pause(0.1)
        if counter == steps - 1:
            plt.pause(0.1)
        else:
            plt.pause(0.1)
            plt.clf()
            
    counter+=1
   
if Animate == False:
    plt.plot(Pos1[:,0],Pos1[:,1], color='purple')
    plt.plot(Pos1[-2:,0],Pos1[-2:,1], color='green')


plt.show()
