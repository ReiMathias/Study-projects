import numpy as np

v0x = 0
v0y = 0
#constants
g = 9.81
a = 0.0065
alpha = 2.5
T0 = 283
K = 4*10**(-5)
Y0 = 10**4
v = 1640
x0 = 0
y0 = 0
angle = 80
h = 0.05
v*np.cos(np.deg2rad(angle))
v*np.sin(np.deg2rad(angle))


#4th order Runge-Kutta algorithm. Takes in current position, time, step-length and function for derivative and returns position at time t+h for position.
def rungekutta(xn, t, h, f, vx0, vy0):
    v0x = vx0
    v0y = vy0
    k1 = f(t,xn,v0x,v0y)
    k2 = f(t + h/2, xn+ (h/2)*k1,v0x,v0y)
    k3 = f(t + h/2, xn+ (h/2)*k2,v0x,v0y)
    k4 = f(t + h, xn+ (h)*k3,v0x,v0y)
    return xn + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    

def rungekuttaspeed(xn, yn, vxn, vyn, t, h, fvx, fvy):
    k1vx = fvx(t,yn, vxn, vyn)
    k1vy = fvy(t,yn, vxn, vyn)
    k1y = vyn
    
    
    k2vx = fvx(t + h/2, yn + h/2*k1y, vxn+ (h/2)*k1vx, vyn+ (h/2)*k1vy)
    k2vy = fvy(t + h/2, yn + h/2*k1y, vxn+ (h/2)*k1vx, vyn+ (h/2)*k1vy)
    k2y = vyn
    
    k3vx = fvx(t + h/2, yn + h/2*k2y, vxn+ (h/2)*k2vx, vyn+ (h/2)*k2vy)
    k3vy = fvy(t + h/2, yn + h/2*k2y, vxn+ (h/2)*k2vx, vyn+ (h/2)*k2vy)
    k3y = vyn
    
    k4vx = fvx(t + h, yn+ h*k3y, vxn + h*k3vx, vyn + h*k3vy)
    k4vy = fvy(t + h, yn+ h*k3y, vxn + h*k3vx, vyn + h*k3vy)
    k4y = vyn
    
    vx = vxn + (h/6)*(k1vx + 2*k2vx + 2*k3vx + k4vx)
    vy = vyn +(h/6)*(k1vy + 2*k2vy + 2*k3vy + k4vy)
    x = xn + h*vx
    y = yn + h*vy

    return [x,y,vx,vy]
    
    

    
#The derivatives of x, y with no drag forces.
def dx(t,xn, v0x,v0y):
    return v0x

def dy(t,yn,v0x,v0y):
    return v0y - g*t
    
#The derivatives of x, y, vx and vy with drag forces assuming adiabatic atmosphere.
def dvxA(t, yn, vx, vy):
    v = np.sqrt(vx**2 + vy**2)
    return( -((1-a*yn/T0)**(alpha)) * K*v*vx )
    
def dvyA(t, yn, vx, vy):
    v = np.sqrt(vx**2 + vy**2)
    print(a,yn,T0,alpha,K,v,vy,g)
    return( -((1-a*yn/T0)**(alpha)) * K*v*vy - g)
    
    
    
#The derivatives of vx and vy assuming isothermal atmosphere.    
def dvxI(t, yn, vx, vy):
    v = np.sqrt(vx**2 + vy**2)
    return -K*np.exp(yn/Y0)*v*vx
    
def dvyI(t, yn, vx, vy):
    v = np.sqrt(vx**2 + vy**2)
    return -K*np.exp(yn/Y0)*v*vy - g    
    
    
    
    
    
    
