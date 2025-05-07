# -*- coding: utf-8 -*- 
""" 
Created on Sun May  4 14:27:29 2025 

@author: bendo 
""" 

import numpy as np 
from scipy.integrate import solve_ivp 
import matplotlib.pyplot as plt 

# defining needed constants 
jmassmodifier = 1 

G = 4*np.pi**2 
M = 1 
Mj = jmassmodifier*(1/1048)*M 
Me = (1/332950)*M 

def round_sig(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags

def earthjuporb(t, coords): 
    xe, xj, ye, yj, vxe, vxj, vye, vyj = coords 
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
  
    output = [dxedt, dxjdt, dyedt, dyjdt, dvxedt, dvxjdt, dvyedt, dvyjdt] 
    return output 

"""  
base integrator section 

""" 
# initial velocity modifier for exageratoin 
ivm=1#.2555 

# starting speed conditions 
EccentricityE = 0.01671123 
Aphelione = 1.016714 
vye0 = np.sqrt(G*M*((1-EccentricityE)/Aphelione))*ivm 
semimajor2 = 5.204267 
AphelionJ = 4.950429 
EccentricityJ = 0.048775 
vxj0 = np.sqrt(G*M*((1-EccentricityJ)/AphelionJ)) 

# starting conditions xe,xj,ye,yj,vxe,vxj,vye,vyj
q0 = [Aphelione, 0, 0, AphelionJ, 0, vxj0, vye0, 0]
# times 
ti = 0 
tf = 100

# integrator error of order (max_step)^9 (for current step size very small)
sol = solve_ivp(earthjuporb, t_span=[ti, tf], y0=q0, method='Radau', max_step=0.01) 

# defining variables after integration 
xesol, xjsol, yesol, yjsol = sol.y[0], sol.y[1], sol.y[2], sol.y[3] 
np.round(xesol,5)
np.round(xjsol,5)
np.round(yesol,5)
np.round(yjsol,5)

""" 
Orbit analysis and seperation 

""" 
# orbit seperations
p=0 
posorb=list() 
for i in range(2,len(yesol)): 

    if yesol[i]>0 and yesol[i-1]<=0: 
        posorb.append(i) 
        p+=1 
  

#difference in steps, found to be constant when tested with different orbits
delta=posorb[1]-posorb[0] 

r=np.sqrt(xesol**2+yesol**2) # position vector magnitude
myfavourite=posorb[0]-delta  # want to get rid of but due to decimal stuff will add an empty orbit could try fix later but this works 
year=1 
R=round_sig(r,7)

#eccentricity collection of instances 
ecoli=list() 
for i in range (posorb[0],len(xesol),delta): 
    #aphelion return 
    aposition=np.argmax(R[myfavourite:i+1]) 
    print('Orbit {} Aphelion value:\t '.format(year),'{}'.format(R[myfavourite+aposition])) 

    #perihelion return 
    periposition=np.argmin(R[myfavourite:i+1]) 
    print('\t\tPerihelion value:','{}\n'.format(R[myfavourite+periposition])) 

    #eccentricity change 
    rarpRatio=R[myfavourite+aposition]/R[myfavourite+periposition] 
    e=(rarpRatio-1)/(rarpRatio+1) 
    ecoli.append(e) 
    myfavourite=myfavourite+delta 

    year=year+1 

#data for plotting first and final orbits 
xorb1b=np.zeros(posorb[0]) 
yorb1b=np.zeros(posorb[0]) 
for i in range(len(xorb1b)+1): 
    xorb1b[i-1]=xesol[i-1] 
    yorb1b[i-1]=yesol[i-1] 
    

xorb2b=np.zeros(delta) 
yorb2b=np.zeros(delta) 
for i in range(delta): 
    xorb2b[i-1]=xesol[posorb[-2]+i-1] 
    yorb2b[i-1]=yesol[posorb[-2]+i-1] 

xorb1=round_sig(xorb1b, 5)
yorb1=round_sig(yorb1b, 5)
xorb2=round_sig(xorb2b, 5)
yorb2=round_sig(yorb2b, 5)

#semi-major evaluation 
rb=np.sqrt(xorb1**2+yorb1**2) 
Rb=round_sig(rb,5)
semib=(Rb[np.argmax(Rb)]+Rb[np.argmin(Rb)])/2 
Roundedsemib=round_sig(semib,5)
print('Initial Semi-Major Axis: {}'.format(Roundedsemib))

rf=np.sqrt(xorb2**2+yorb2**2) 
Rf=round_sig(rf,5)
semif=(Rf[np.argmax(Rf)]+Rf[np.argmin(Rf)])/2 
print(Rf[np.argmax(Rf)],Rf[np.argmin(Rf)]) 
Roundedsemif=round_sig(semif,5)
print('\nFinal Semi-Major Axis: {}'.format(round_sig(semif,5)))

print('\nInitial Perihelion Position: ', xorb1[[np.argmin(Rb)]],yorb1[[np.argmin(Rb)]]) 
print('Final Perihelion Position: ', xorb2[[np.argmin(Rf)]],yorb2[[np.argmin(Rf)]]) 

k=round_sig(ecoli[0],5)
l=round_sig(ecoli[-1],5)

print('\nInitial Eccentricity: {}'.format(k)) 
print('Final Eccentricity: {}'.format(l)) 

print('Delta Semi-major: {}'.format(Roundedsemif-Roundedsemib)) 
print('Delta Eccentricity: {}'.format(round_sig(l-k, 5))) 

""" 
plotting 

""" 

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
plt.title("Earth Jupiter Orbit ivp, {}mj, {}ve, {}years".format(jmassmodifier,ivm,tf)) 
plt.show() 

plt.subplot(1,1,1) 
plt.plot(xorb1,yorb1, label="earthorb initial") 
plt.plot(xorb2,yorb2, label="earthorb final", color='r') 
plt.xlabel("x") 
plt.ylabel("y") 
plt.xlim(0, 1.2) 
plt.ylim(0, 1) 
plt.legend() 
plt.title("Earth close up, Orbit ivp,{}mj, {}ve,{}years".format(jmassmodifier,ivm,tf)) 
plt.show() 

plt.plot(0,0,'ro', label='sun') 
plt.plot(xorb1,yorb1, label="earthorb initial") 
plt.plot(xorb2,yorb2, label="earthorb final", color='r') 
plt.grid(visible=True, which='major', color='#DDDDDD', linestyle='-') 
plt.grid(visible=True, which='minor', color='#EEEEEE', linestyle='--') 
plt.minorticks_on() 
plt.xlabel("x") 
plt.ylabel("y") 
plt.legend() 
plt.title("Earth initial and final comp, {}mj, {}ve, {}years".format(jmassmodifier,ivm,tf)) 
plt.show() 

plt.plot(range(ti,tf), ecoli) 
plt.grid(visible=True, which='major', color='#DDDDDD', linestyle='-') 
plt.grid(visible=True, which='minor', color='#EEEEEE', linestyle='--') 
plt.minorticks_on() 
plt.xlabel("years") 
plt.ylabel("eccentricity") 
plt.title("Earth eccentricity, {}mj, {}ve {}years".format(jmassmodifier,ivm,tf)) 
plt.show() 


""" 
error was determined based off order of integration method being to 6 decimal
places and then worst sig figs from integration values of 5 sig fig
this shows constant semimajor as expected

""" 

 