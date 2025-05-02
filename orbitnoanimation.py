# -*- coding: utf-8 -*-
"""
Created on Thu May  1 15:47:23 2025

@author: bendo
"""


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation # see bottom


def earthjuporb(t, coords):
    xe, xj, ye, yj, vxe, vxj, vye, vyj = coords

    # defining needed constants
    G = 4*np.pi**2
    M = 1
    Mj = 10*(1/1048)*M
    Me = (1/332950)*M

    re = np.sqrt(xe**2+ye**2)
    rj = np.sqrt(xj**2+yj**2)
    rej = np.sqrt((xe-xj)**2+(ye-yj)**2)

    # speed derivatives
    dxedt = vxe
    dxjdt = vxj
    dyedt = vye
    dyjdt = vyj

    # speed derivative forms of acceleration
    dvxedt = -G*(M*xe/((re)**3)+Mj*(xj-xe)/((rej)**3))
    dvyedt = -G*(M*ye/((re)**3)+Mj*(yj-ye)/((rej)**3))
    dvxjdt = -G*(M*xj/((rj)**3)-Me*(xj-xe)/((rej)**3))
    dvyjdt = -G*(M*yj/((rj)**3)-Me*(yj-ye)/((rej)**3))

    # derivatives
    output = [dxedt, dxjdt, dyedt, dyjdt, dvxedt, dvxjdt, dvyedt, dvyjdt]
    return output


# starting speed conditions
G = 4*np.pi**2 #better than previous decimal approximation, which was causing the orbits to be too fast and add extra
Eccentricity = 0.01671123
M = 1
Aphelion = 1.016714
vye0 = np.sqrt(G*M*(2/Aphelion-1/1))*1.2555
semimajor = 5.204267
Aphelion2 = 4.950429
Eccentricity2 = 0.048775
vxj0 = np.sqrt(G*M*(2/Aphelion2-1/semimajor))

# starting conditions xe,xj,ye,yj,vxe,vxj,vye,vyj
start = [Aphelion, 0, 0, Aphelion2, 0, vxj0, vye0, 0]

# time array
t = np.linspace(0, 110.859, 550)

# integrator, should remove time array leftover from odeint and euler remnants
sol = solve_ivp(earthjuporb, t_span=[t[0], t[len(
    t)-1]], y0=[Aphelion, 0, 0, Aphelion2, 0, vxj0, vye0, 0], method='Radau', max_step=0.01)

# sanity check to see inner workings can be hidden away when happy
print(sol.y[0])

# defining variables after integration
xesol, xjsol, yesol, yjsol = sol.y[0], sol.y[1], sol.y[2], sol.y[3]







#sanity check definitions, Im commenting a lil backwards
p=0
posorb=list()
for i in range(2,len(yesol)):
    if yesol[i]>=0 and yesol[i-1]<0:
        print('\n',i)
        print(yesol[i],'\n')
        posorb.append(i)
        p+=1
    
#sanity checks to be removed upon completion but show things are working
print(p)
print(posorb,'\n\n')

#difference in steps, found to be constant which is nice
delta=posorb[1]-posorb[0]

# defining the actual position magnitude from focus
# have tried getting rid of q but for some reason it adds a fourth orbit?? not cool
# will fix in clean up but want this to actually work
R=np.sqrt(xesol**2+yesol**2)
myfavourite=posorb[0]-delta  # want to get rid of but due to decimal stuff will add an empty orbit could try fix later but this works
year=1
for i in range (posorb[0],len(xesol),delta):
    
    #aphelion return
    aposition=np.argmax(R[myfavourite:i+1])
    print('Orbit {} Aphelion value:\t '.format(year),R[myfavourite+aposition])
    
    #perihelion return
    periposition=np.argmin(R[myfavourite:i+1])
    print('\t\tPerihelion value:',R[myfavourite+periposition],'\n')
    myfavourite=myfavourite+delta
    year=year+1




xorb1=np.zeros(posorb[0])
yorb1=np.zeros(posorb[0])
for i in range(len(xorb1)+1):
    xorb1[i-1]=xesol[i-1]
    yorb1[i-1]=yesol[i-1]

xorb2=np.zeros(delta)
yorb2=np.zeros(delta)



for i in range(posorb[-1]-posorb[-2]+1):
    xorb2[i-1]=xesol[posorb[-2]+i-1]
    yorb2[i-1]=yesol[posorb[-2]+i-1]


jmassmodifier = 10

plt.plot(0,0,'ro', label='sun')
plt.plot(xesol,yesol, label="earthorb")
plt.plot(xjsol,yjsol, label="juporb")
#fig.patch.set_facecolor('black')
plt.grid(visible=True, which='major', color='#DDDDDD', linestyle='-')
plt.grid(visible=True, which='minor', color='#EEEEEE', linestyle='--')
plt.minorticks_on()
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.legend()
plt.title("Earth Jupiter Orbit ivp, {}mj".format(jmassmodifier))
plt.show()


plt.subplot(1,1,1)
plt.plot(xorb1,yorb1, label="earthorb initial")
plt.plot(xorb2,yorb2, label="earthorb final", color='r')
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.title("Earth close up, Orbit ivp, {}mj".format(jmassmodifier))
plt.show()


plt.plot(0,0,'ro', label='sun')
plt.plot(xorb1,yorb1, label="earthorb initial")
plt.plot(xorb2,yorb2, label="earthorb final", color='r')
#fig.patch.set_facecolor('black')
plt.grid(visible=True, which='major', color='#DDDDDD', linestyle='-')
plt.grid(visible=True, which='minor', color='#EEEEEE', linestyle='--')
plt.minorticks_on()
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.title("Earth initial and final comp, {}mj".format(jmassmodifier))
plt.show()


"""
how to make this an animation for later, rather than doing stupid annoying ball
flying around just show orbits changing with each orbit being a new frame should
be real easy and more visually useful than the ball
""" 