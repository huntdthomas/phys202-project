import numpy as np
from scipy.integrate import odeint

Direct = 1
parabolic = 1

if parabolic != 1:
    epsilon = 0.2
    
gamma = 4.49933e4 #units of (kpc^3) / (M_sun^10)(years^9)^2    Derived by James Amarel

# M and S are the masses of the central galaxy and disrupting galaxy
if parabolic == 1:
    M = 1
    S = 1
else:
    M = 5
    S = 10

# atol and rtol are the error tolerances for the odeint solver.
atol,rtol = 1e-6,1e-6

#Disrupting Galaxy S initial position in y
if parabolic == 1:
    Sy = 40
else:
    Sy = 30*(1+epsilon)
# tmax is the maximum time the simulation will run to.
# t is a numpy array with equally spaced values of time starting from 0 and enting at tmax.

if parabolic == 1:
    tmax = 1.0
else:
    tmax = 5.0
t = np.linspace(0,tmax,1000)


def ShellFill(number,radius,Direct):
    ic = []
    thetastep = 2*np.pi/number
    theta = np.arange(0,2*np.pi,thetastep)
    s = [np.sin(k) for k in theta]
    c = [np.cos(k) for k in theta]
    if Direct == 1:
        v = np.sqrt(gamma*M/radius)
    else:
        v = -np.sqrt(gamma*M/radius)
    for i in range(len(theta)):
        newic = [radius*s[i],v*c[i],radius*c[i],-v*s[i]]
        ic.append(newic)
    return ic

def makeN(Direct):  
    """
    r0: A base radius to start the stars at.
    thetastep, theta: These are the step size, defined by number of stars, and angular   position
                    of the stars. Note how the stars are set up to have equally spaced angles.
    ic: The array of initial conditions for the positions and velocities of the stars:
            of the form [[x1,dx1,y1,dy1],[x2,dx2,y2,dy2]...[xf,dxf,yf,dyf]] for N*shells stars.

    The paper has radii [0.2,0.3,0.4,0.5,0.6]*rmin, with [12,18,24,30,36] stars at each level
    """
    r0 = 25
    N = np.array([12,18,24,30,36])
    shells = np.array([0.20,0.30,0.40,0.50,0.60])*r0
    ic = []
    
    I = list(zip(N,shells))
    for i in I:
        newic = ShellFill(i[0],i[1],Direct)
        for i in newic:
            ic.append(i)
    return ic




def parab(y):
    rmin = 25
    x = rmin - (y**2)/(4*rmin)
    R = np.sqrt(x**2 + y**2)
    v = np.sqrt(2*gamma*(M+S)/R)
    angle = np.arctan(50/y)
    if y == 0:
        vx,vy = 0,-v
    else:
        vy = -v*np.sin(angle)
        vx = v*np.cos(angle)
    return [x,vx,y,vy]




def ellip(y,epsilon):
    
    rmin = 30
    c = rmin*(1+epsilon)
    a = c/(1-epsilon**2)
    b = c/np.sqrt(1 - epsilon**2)
    d = a*epsilon
    x = d - a*np.sqrt(-(y**2)/(b**2) + 1)
    R = np.sqrt(x**2 + y**2)
    v = np.sqrt(gamma*(M+S)*((2/R)-(1/a)))
    angle = np.arctan(-(1/a*y)*np.sqrt((-(y**2)/(b**2)) + 1))
    
    vy = -v*np.sin(angle)
    vx = v*np.cos(angle)

    return [x,vx,y,vy]




def Galaxy(parabolic,Direct,Sy):
    """
    X,Y: The relative X and Y positions of the disrupting galaxy to the central galaxy.
    R: The relative distance of the disrupting galaxy to the central galaxy.
    Vs: The velocity of the disrupting galaxy
    Sinfo: The array of positions, velocities of the disrupting galaxy: [x,dx,y,dy]
    """
    if parabolic == 1:
        Sinfo = parab(Sy)
    else:
        Sinfo = ellip(Sy,epsilon)
    # calls makeN() to pull the initial conditions for the stars to ic
    ic = makeN(Direct)
    # this for loop pulls open each array in ic (for each individual star) and appends to it the
    # initial conditions for the galaxy S, Sinfo.
    for i in ic:
        for j in range(len(Sinfo)):
            i.append(Sinfo[j])
    ic = np.array(ic)
    return ic




def derivs(rvec,t):
    """
    rvec = [x,dx,y,dy,X,dX,Y,dY]
    
    x:  x component of relative position of single star to central galaxy.
    dx: x component of relative velocity of single star to central galaxy.
    y:  y component of relative position of single star to central galaxy.
    dy: y component of relative velocity of single star to central galaxy.
    
    X:  x component of relative position of disrupting galaxy to central galaxy.
    dX: x component of relative velocity of disrupting galaxy to central galaxy.
    Y:  y component of relative position of disrupting galaxy to central galaxy.
    dY: y component of relative velocity of disrupting galaxy to central galaxy.
    """
    
    x,y,X,Y = rvec[0],rvec[2],rvec[4],rvec[6]   
    dx,dy,dX,dY = rvec[1],rvec[3],rvec[5],rvec[7]
    
    """
    r: The radial distance from the single star to the central galaxy.
    R: The radial distance from the central galaxy to the disrupting galaxy.
    """
    r = np.sqrt(x**2 + y**2)
    R = np.sqrt(X**2 + Y**2)
    
    """
    ddx,ddy: The x and y components of the radial acceleration of the single star, respectively.
    ddX,ddY: The X and Y components of the radial acceleration of the disrupting galaxy, respectively.
    """
    
    ddx = -gamma*((M/(r**3))*x - (S/(np.sqrt((X-x)**2 + (Y-y)**2))**3)*(X-x) + (S/(R**3)*X))
    ddy = -gamma*((M/(r**3))*y - (S/(np.sqrt((X-x)**2 + (Y-y)**2))**3)*(Y-y) + (S/(R**3)*Y))
    
    ddX = -gamma*((M+S)/(R**3))*X
    ddY = -gamma*((M+S)/(R**3))*Y
                  
    return np.array([dx,ddx,dy,ddy,dX,ddX,dY,ddY])



def get_solution(parabolic,Direct,Sy):
    ic = Galaxy(parabolic,Direct,Sy)
    #these are the x and y positions/velocities for all the stars
    x,dx,y,dy = [],[],[],[]
    for k in range(len(ic)):
        soln = odeint(derivs,ic[k],t,atol=atol,rtol=rtol)
        x.append(np.array([j[0] for j in soln]))
        y.append(np.array([j[2] for j in soln]))     
        dx.append(np.array([j[1] for j in soln]))
        dy.append(np.array([j[3] for j in soln]))
    
    # note how the initial conditions for the disrupting galaxy can be defined for any iteration
    # as you iterate through k, soln is left with the last value, but the Sinfo components are the same
    #the following are the X and Y position/velocities of the disrupting galaxy
    X = np.array([j[4] for j in soln])
    Y = np.array([j[6] for j in soln])  
    dX = np.array([j[5] for j in soln])
    dY = np.array([j[7] for j in soln])
    
    return x,dx,y,dy,X,dX,Y,dY


def energy(parabolic,Direct,Sy):
    x,dx,y,dy,X,dX,Y,dY = get_solution(parabolic,Direct,Sy)
    
    kinetic = 0.5*S*(dX**2+dY**2)  
    potential = -gamma*M*S/np.sqrt((X**2+Y**2))
    E = (kinetic + potential)
    
    return E

#assign solutions to a giant array. Pass this to plot_motion(see below).
def Save(parabolic,Direct,Sy):
    soln = get_solution(parabolic,Direct,Sy)
    if parabolic == 1:
        if Direct == 1:
            np.savez("ParabDirect", x= soln[0],dx = soln[1], y = soln[2],dy = soln[3],X = soln[4],dX = soln[5],Y = soln[6],dY = soln[7])
        else:
            np.savez("ParabRetro", x= soln[0],dx = soln[1], y = soln[2],dy = soln[3],X = soln[4],dX = soln[5],Y = soln[6],dY = soln[7])
    else:
        if Direct == 1:
            np.savez("EllipDirect", x= soln[0],dx = soln[1], y = soln[2],dy = soln[3],X = soln[4],dX = soln[5],Y = soln[6],dY = soln[7])
        else:
            np.savez("EllipRetro", x= soln[0],dx = soln[1], y = soln[2],dy = soln[3],X = soln[4],dX = soln[5],Y = soln[6],dY = soln[7])