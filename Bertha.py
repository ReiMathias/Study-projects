import functions as fun
import numpy as np
from matplotlib import pyplot as plt

g = 9.81
#Initialization
v = 1640
x0 = 0
y0 = 0
angle = 80
h = 0.05
v0x = v*np.cos(np.deg2rad(angle))
v0y = v*np.sin(np.deg2rad(angle))

def yesDragAdiabat():
    t = 0
    x = np.array([x0])
    y = np.array([y0])
    vx = np.array([v0x])
    vy = np.array([v0y])
    yi = y0
    while yi>= 0:
        t+=0
        rk = fun.rungekuttaspeed(x[-1], y[-1], vx[-1], vy[-1],t, h, fun.dvxA, fun.dvyA)
        yi = rk[1]
        x = np.append(x, rk[0])
        y = np.append(y, rk[1])
        vx = np.append(vx, rk[2])
        vy = np.append(vy, rk[3])
    return [x,y,x[-1], np.amax(y)]



list = yesDragAdiabat()
print('Range of Big Bertha: ', list[2], ' Height of Big Bertha: ', list[3])

plt.plot(list[0],list[1], label = "Projectile of Big Bertha")

plt.xlabel("Length [m]")
plt.ylabel("Height [m]")
plt.legend(loc = 1)
plt.grid()
plt.show()