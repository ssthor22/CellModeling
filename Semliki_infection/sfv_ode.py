# Author: Seth Thor
# Zan lab rotation 
# Date created: Nov 2014

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

#Constants
tstart = 0
tend = 7200 #2 hours

#Rates
kaC = 0.00035	#(/s) Note: ka = 5.2E-9 mL/cell/min
ke = 0.00065	#(/s) endocytosis
kfusn = 0.0023	#(/s) fusion into endosomal membrane
ktran = kfusn	#(/s) transport inactive viruses to lysosomes, equal when 

#Initial species count
Vex = 1000
Vi = 0
Vs = 100
V_endotot = 0
Vcyt = 0
V_o = [Vex, Vi, Vs, V_endotot, Vcyt]

#ODE
def dVdt(s, t):
    Vex_i = s[0]
    Vi_i = s[1]
    Vs_i = s[2]
    V_endotot_i = s[3]
    Vcyt_i = s[4]
    #rate eqs
    dVex_dt = -kaC*Vex_i
    dVi_dt = ke*Vs_i
    dVs_dt = kaC*Vex_i - ke*Vs_i
    dV_endotot_dt = ke*Vs_i - kfusn*V_endotot_i - ktran*V_endotot_i
    dVcyt_dt = kfusn*V_endotot_i
    return [dVex_dt, dVi_dt, dVs_dt, dV_endotot_dt, dVcyt_dt]
    
#Solve
t = np.linspace(tstart, tend, 100000)
soln = spi.odeint(dVdt, V_o, t)

#Plot
plt.figure()
plt.plot(t, soln[:,0], label="Vex")
plt.plot(t, soln[:,1], label="Vi")
plt.plot(t, soln[:,2], label="Vs")
plt.plot(t, soln[:,3], label="V_endotot")
plt.plot(t, soln[:,4], label="Vcyt")
plt.xlabel("Time (s)")
plt.ylabel("Virus count")
plt.title("Pre-bound + Extracellular Viruses")
plt.legend()
plt.savefig("sfv_ode2.png")

print "Done"