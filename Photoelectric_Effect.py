# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 18:48:51 2024

@author: Michael Younger
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import h, k, c 

el_charge=1.602e-19
sr_on_earth = 6.794e-5
sun_temp = 5775

Ens = np.linspace(43, 4429, num=4386) / 1e3
lambdas = h*c / (Ens * el_charge) * 1e9

def planck_e(energy, T):
    
    energy=energy*el_charge  # convert to Joules
    a2 = 2.0/(h*c)**2 * (el_charge/h)
    b2 = energy/(k*T)
    
    Intensity2 = a2*energy**3 /( (np.e**b2 - 1.0) )
    
    return Intensity2 #back in units of (W/(sr*m^2*eV)

def planck_l(wav, T):
    
    wav=wav*1e-9  # convert to meters
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    
    Intensity = a/( (wav**5) * (np.e**b - 1.0) )
    
    return Intensity #back in units of W/(sr*m^3)

#calculate what to plot
AM0_l = planck_l(lambdas, sun_temp) * sr_on_earth * 1e-9
AM0_e = planck_e(Ens, sun_temp) * sr_on_earth
AM0_e_nph = AM0_e / Ens
AM0_l_nph = AM0_l / Ens

# Plot stuff
plt.plot(lambdas[300:], AM0_l[300:], label=f'T={sun_temp}',  color='red')  #The first 300 cause the plot to skew a lot, creates better looking plot

# Annotate plots
plt.xlabel(r'$\lambda$ (nm)')
plt.ylabel('Spectral Irradiance (W/(m^2*nm))')
#plt.legend()
plt.show()

plt.plot(Ens, AM0_e, color='red')

# Annotate plots
plt.xlabel('Energy (eV)')
plt.ylabel('Spectral Energy Density (W/(m^2*eV))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

plt.plot(lambdas[300:], AM0_l_nph[300:], color='red')

# Annotate plots
plt.xlabel(r'$\lambda$ (nm)')
plt.ylabel('Photon Flux (Photons/(s*m^2*nm))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

plt.plot(Ens, AM0_e_nph, color='red')

# Annotate plots
plt.xlabel('Energy (eV)')
plt.ylabel('Photon Flux (Photons/(s*m^2*eV))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

#
# Calculate and compare W/m^2
#

print('The integrated spectral intensity vs wavelength for T={:d} is\n {:.2f} W/m²' 
      .format(sun_temp, np.trapz(AM0_l, x = -1 * lambdas)))

print('The integrated spectral intensity vs energy for T={:d} is\n {:.2f} W/m²' 
      .format(sun_temp, np.trapz(AM0_e, x=Ens)))

print('The integrated spectral intensity based on photon flux vs wavelength for T={:d} is\n {:.2f} W/m²'
      .format(sun_temp, np.trapz(AM0_l_nph * Ens, x = -1 * lambdas)))

print('The integrated spectral intensity based on photon flux vs energy for T={:d} is\n {:.2f} W/m²' 
      .format(sun_temp, np.trapz(AM0_e_nph * Ens, x=Ens)))

#
# Determine power harvested by a Silicon solar cell
#

# find the closest eV value in the array
# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_indices(array, values):
    array = np.asarray(array)
    indices = []
    for value in values:
        idx = (np.abs(array - value)).argmin()
        indices.append(idx)
    return indices
#
#  Find array index for which E = 1100 meV
#
array_index_Si = find_nearest_indices(Ens, values=[1.1])[0]

#
# Integrate over energy array, starting at Egap(Si), multiplying with Egap instead of photon energy
#
received_power_density = np.trapz(AM0_e_nph * Ens, x=Ens)
harvested_power_density_Si = np.trapz(AM0_e_nph[array_index_Si:]*(1.1), x=Ens[array_index_Si:])


print('The power harvested from a Si solar cell is {:.2f} W/m², corresponding to an efficiency of {:.2%} '
      .format(harvested_power_density_Si,harvested_power_density_Si/received_power_density))


Vg = np.linspace(0.5, 2.5, 100)
efficiency = []
array_indeces = np.array(find_nearest_indices(Ens, Vg))

for n, idx in enumerate(array_indeces):
    efficiency.append(np.trapz(AM0_e_nph[idx:]*(Vg[n]), x=Ens[idx:]) / np.trapz(AM0_e_nph * Ens, x=Ens))



plt.plot(Vg, np.array(efficiency), color='black')

# Annotate plots
plt.xlabel('Band gap energy')
plt.ylabel('Efficiency')
#plt.xscale("log")
#plt.yscale("log")
plt.show()



