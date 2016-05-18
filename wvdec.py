# Copyright (C) 2014,2015,2016 Joern Callies
#
# This file is part of WV-dec.
#
# WV-dec is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WV-dec is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WV-dec. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

# Two examples for loading in spectra are given below. The variables are:
# k -  wavenumber
# Su - along-track velocity spectrum
# Sv - across-track velocity spectrum
# Sw - verical velocity spectrum
# Sb - potential energy spectrum (buoyancy spectrum normalized by N^2)

# This reads in a set of spectra saved in a numpy archive.
f = np.load('start08_strat.npz')
k = f['k']
Su = f['Su']
Sv = f['Sv']
Sw = f['Sw']
Sb = f['Sb']

# This reads in a set of spectra saved in a MATLAB file. [uncomment for use]
#f = sio.loadmat('start08_strat.mat')
#k = f['k'][0]
#Su = f['Su'][0]
#Sv = f['Sv'][0]
#Sw = f['Sw'][0]
#Sb = f['Sb'][0]

# Helmholtz decomposition
# This decomposes the horizontal kinetic energy K = 1/2(Su+Sv) into a rotational
# part Kpsi and a divergent part Kphi. This formulation uses equations (2.11)
# and (2.12) of Lindborg (2015, JFM), which follow directly from (2.27), (2.30),
# and (2.31) of Buehler et al. (2014, JFM). The integration is performed using a
# logarithmic integration variable.
K = (Su+Sv)/2
Kphi = np.empty(k.size)
Kpsi = np.empty(k.size)
for i in range(k.size):
    Kphi[i] = (Su[i]-np.trapz(k[i:]*(Sv[i:]-Su[i:]), x=np.log(k[i:]))/k[i])/2
    Kpsi[i] = (Sv[i]+np.trapz(k[i:]*(Sv[i:]-Su[i:]), x=np.log(k[i:]))/k[i])/2

# Hydrostatic wave-vortex decomposition
# This decomposes the total energy E = 1/2(Su+Sv+Sb) into a wave component Ew
# and a residual vortical component Ev assuming hydrostatic waves. [uncomment
# for use]
#E = (Su+Sv+Sb)/2
#Ew = 2*Kphi
#Ev = E-Ew

# Nonhydrostatic wave-vortex decomposition
# This decomposes the total energy E = 1/2(Su+Sv+Sw+Sb) into a wave component Ew
# and a residual vortical component Ev assuming hydrostatic waves.
E = (Su+Sv+Sw+Sb)/2
Ew = 2*Kphi+Sw
Ev = E-Ew

# Plot the Helmholtz decomposition result.
plt.figure()
plt.loglog(k, K)
plt.loglog(k, Kpsi)
plt.loglog(k, Kphi)
plt.xlabel('wavenumber [rad/m]')
plt.ylabel('power spectral density [kg/s^2]')
plt.legend(['horizontal kinetic energy', 'rotational component', 'divergent component'], loc=3)

# Plot the wave-vortex decomposition result.
plt.figure()
plt.loglog(k, E)
plt.loglog(k, Ew)
plt.loglog(k, Ev)
plt.xlabel('wavenumber [rad/m]')
plt.ylabel('power spectral density [kg/s^2]')
plt.legend(['total energy', 'diagnosed wave component', 'residual vortical component'], loc=3)

plt.show()
