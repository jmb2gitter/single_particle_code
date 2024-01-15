""" TPC spatial profiles estimator. The magnetic configuration is of the one
    magnetic bottle, for which the magnetic intensity at its center is equal
    to the equatorial magnetic field of the Lshell especified, and the magnetic
    intensity on its end points is equal to the value of the field at the
    intesection with the earths surface. """


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

''' ----- initial parameters and physical constants ----- '''
qe = -1.602176462e-19   # electron charge (C)
me = 9.10938188e-31     # electron mass (kg)
Be = 3.11e-5            # equatorial magnetic field (T)
Re = 6.371e6            # earth radius (km)
L = 4                   # magnetic field line

''' ----- simulation_environment_setup ----- '''

# latitude at which the magnetic field line intersects the earths surface (rad)
li = math.acos(1 / math.sqrt(L))
# estimation of the magnetic field line length (Re)
si = 0.5*L*Re*(math.sin(li)*math.sqrt(1 + 3*math.sin(li)**2) +
               math.log(math.sqrt(3)*math.sin(li) +
                        math.sqrt(1 + 3*math.sin(li)**2))/math.sqrt(3))
# spatial discretization units
ns = 500
# spatial discretization interval
ds = 2*si/ns
# spatial domain declaration & initialization (Re)
s = np.linspace(-si, si, 501)
# time discretization interval
dt = .01

# magnetic field setup - magnetic bottle B = Beq (1 + (s/so)^2)
# magnetic field intensity at que equator:     Beq = Be/L^3
Beq = Be/L**3
# magnetic field intensity at the ionosphere:  Bi = Beq * sqrt(1 + 3*sin(li)^2)/cos(li)^6
Bi =  Beq * math.sqrt(1 + 3*math.sin(li)**2)/math.cos(li)**6
# so = si / sqrt(Bi/Beq - 1)
so = si/math.sqrt(Bi/Beq - 1)
# magnetic intensity declaration & initialization (T):
magneticField = np.array([Beq*(1 + (x/so)**2) for x in s])
# magnetic field profile from simulation code
spatial_profiles = pd.read_csv('spatial_profiles.csv')

plt.plot(spatial_profiles['position (Re)']/1e3, spatial_profiles['magnetic field (T)'], 'b-', label='simulation')
plt.plot(s/1e3, magneticField, 'g-', label=r"$B(s) = B_0\left (1 + \frac{s^2}{s_0^2} \right)$")
plt.legend()
plt.xlabel('distance (km)')
plt.ylabel('magnetic intensity (T)')
#plasma_initialization()
#wave_initialization()

"""for i in range(4):
    for p in range(15, 90, 15):
        
un_electron = pd.read_csv("electron20.csv")
plt.plot(un_electron[' v_perp(m/s)'], un_electron[' v_par(m/s)'], '.')
"""
electron = pd.read_csv("electron20.csv")
electron['time(s)'] = dt * electron.index
electron.set_index('time(s)', inplace=True)
electron["energy(eV)"] = .5 * me * (electron[' v_par(m/s)']**2 + electron[' v_perp(m/s)']**2 ) / (-qe)

fig, (ax1, ax2) = plt.subplots(2)
electron.plot('s(m)', ' v_par(m/s)', ax=ax1)
ax1.set_ylabel("$v_\parallel$ (m/s)")
ax1.set_xlabel("$s$ (m)")
electron.plot(y='s(m)', ax=ax2)
electron.plot(y='energy(eV)', ax=ax2, secondary_y=True)
plt.ylim((990, 1010))