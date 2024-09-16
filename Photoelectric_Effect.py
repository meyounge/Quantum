# -*- coding: utf-8 -*-
"""
Created on Aug 20 11:56:31 2024

@author: mgoryll

Copyright: Arizona State University

Python script to visualize the quantum mechanical nature of the 
Thermodynamic Efficiency Limit for Photovoltaic Devices

"""

#
# Plot Planck's energy distribution vs wavelength for different temperatures
#
# Source: https://web.archive.org/web/20221105082528/https://dpotoyan.github.io/Chem324/python-intro.html
#

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import h, k, c 

el_charge=1.602e-19
sr_on_earth = 6.794e-5

# Define range of wavelengths (nanometers)
#lambdas = np.linspace(280, 4e3, num=2000) # NREL spacing
lambdas = np.linspace(280, 3e3, num=2000)
sun_temp = 5775

# Define convenient function for computing spectral intensity as a function of T
def planck_l(wav, T):
    
    wav=wav*1e-9  # convert to meters
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    
    Intensity = a/( (wav**5) * (np.e**b - 1.0) )
    
    return Intensity #back in units of W/(sr*m^3)

# Plot intensities, vary T
plt.plot(lambdas, planck_l(lambdas,sun_temp) * sr_on_earth * 1e-9, label=f'T={sun_temp}',  color='red')  

# Annotate plots
plt.xlabel(r'$\lambda$ (nm)')
plt.ylabel('Spectral Irradiance (W/(m^2*nm))')
#plt.legend()
plt.show()

#
# Determine the received intensity on Earth (ignoring atmospheric effects)
# 
AM0_l = planck_l(lambdas, sun_temp) * sr_on_earth * 1e-9 # eliminate sr and 1/nm

#
# Calculate the integral over wavelengths
#
print('The integrated spectral intensity for T={:d} is {:.2f} W/m²' .format(sun_temp, np.trapz(AM0_l)))

#
# Plot Planck's intensity distribution vs frequency for different temperatures
#
'''
# Define range of frequencies (THz)
#nus = np.linspace(75, 1.07e3, num=2000) # NREL spacing
nus = np.linspace(100, 1.07e3, num=970)

# Define convenient function for computing spectral energy density as a function of T
def planck_f(freq, T):
    
    freq=freq*1e12  # convert to THz
    a1 = 2.0*h/c**2
    b1 = h*freq/(k*T)
    
    Intensity1 = a1*freq**3 /( (np.e**b1 - 1.0) )
    
    return Intensity1 # back in units of (W/(sr*m^2*Hz)

# Plot intensities, vary T
plt.plot(nus, planck_f(nus,sun_temp), color='red')  

# Annotate plots
plt.xlabel('Frequency (THz)')
plt.ylabel('Spectral Energy Density (W/(sr*m^2*Hz))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

#
# Determine the received intensity on Earth (ignoring atmospheric effects)
# 
AM0_f = planck_f(nus,sun_temp)*6.794e-5*1e12 # eliminate sr and 1/THz

#
# Calculate the integral over frequencies
#
print('The integrated spectral intensity for T={:d} is {:.2f} W/m²' .format(sun_temp, np.trapz(AM0_f)))

'''
#
# Plot Planck's intensity distribution vs energy for different temperatures
#

# Define range of energies (electron volts)
Ens = np.linspace(43, 4429, num=4386) / 1e3
#lambdas = h*c / (Ens * el_charge) * 1e9 
#Ens = h*c / (lambdas / 1e9) / el_charge

# Define convenient function for computing spectral energy density as a function of T
# Energy(eV) = h*freq/el_charge <=> freq = Energy*el_charge/h <=> df = dE * el_charge/h
def planck_e(energy, T):
    
    energy=energy*el_charge  # convert to Joules
    a2 = 2.0/(h*c)**2 * (el_charge/h)
    b2 = energy/(k*T)
    
    Intensity2 = a2*energy**3 /( (np.e**b2 - 1.0) )
    
    return Intensity2 #back in units of (W/(sr*m^2*eV)

# Plot intensities, vary T
plt.plot(Ens, planck_e(Ens,sun_temp) * sr_on_earth, color='red')

# Annotate plots
plt.xlabel('Energy (eV)')
plt.ylabel('Spectral Energy Density (W/(m^2*eV))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

#
# Determine the received intensity on Earth (ignoring atmospheric effects)
# 
#AM0_e = planck_e(Ens,6000)*6.794e-5*el_charge/(1e3*h) # eliminate sr and 1/eV
AM0_e = planck_e(Ens, sun_temp) * sr_on_earth # eliminate sr

#
# Calculate the integral over energies
#
print('The integrated spectral intensity for T={:d} is {:.2f} W/m²' .format(sun_temp, np.trapz(AM0_e) / 1e3))

#
# Convert the intensity array to a photon flux array
#
# Divide by photon energy in Joules 
# Re-normalize bin from meV to eV by multiplying by 1e3
AM0_e_nph = AM0_e / Ens /el_charge / 1e3

plt.plot(Ens, AM0_e / Ens, color='red')

# Annotate plots
plt.xlabel('Energy (eV)')
plt.ylabel('Photon Flux (Photons/(s*m^2*eV))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

plt.plot(lambdas, planck_l(lambdas,sun_temp) * sr_on_earth * 1e-9 / Ens, color='red')

# Annotate plots
plt.xlabel(r'$\lambda$ (nm)')
plt.ylabel('Photon Flux (Photons/(s*m^2*nm))')
#plt.xscale("log")
#plt.yscale("log")
plt.show()

#
# Calculate the integral over energies based on photon flux
#
print('The integrated spectral intensity based on photon flux for T={:d} is {:.2f} W/m²' 
      .format(sun_temp, np.trapz(AM0_e_nph*(Ens*el_charge))))
received_power_density = np.trapz(AM0_e_nph*(Ens*el_charge))

#
# Determine power harvested by a Silicon solar cell
#

# find the closest eV value in the array
# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#
#  Find array index for which E = 1100 meV
#
array_index_Si = find_nearest_idx(Ens, value=1.1)

#
# Integrate over energy array, starting at Egap(Si), multiplying with Egap instead of photon energy
#
harvested_power_density_Si = np.trapz(AM0_e_nph[array_index_Si:]*(1.1*el_charge))

print('The power harvested from a Si solar cell is {:.2f} W/m², corresponding to an efficiency of {:.2f} '
      .format(harvested_power_density_Si,harvested_power_density_Si/received_power_density))

