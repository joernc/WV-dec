% Copyright (C) 2014,2015,2016 Joern Callies
%
% This file is part of WV-dec.
%
% WV-dec is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% WV-dec is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with WV-dec. If not, see <http://www.gnu.org/licenses/>.

% This reads in a set of spectra. The variables are:
% k -  wavenumber
% Su - along-track velocity spectrum
% Sv - across-track velocity spectrum
% Sw - verical velocity spectrum
% Sb - potential energy spectrum (buoyancy spectrum normalized by N^2)
load start08_strat.mat

% Helmholtz decomposition
% This decomposes the horizontal kinetic energy K = 1/2(Su+Sv) into a rotational
% part Kpsi and a divergent part Kphi. The integration is performed using a
% logarithmic integration variable.
K = (Su+Sv)/2;
Kphi = zeros(size(k));
Kpsi = zeros(size(k));
for i = 1:length(k)
  Kphi(i) = (Su(i)-trapz(log(k(i:end)), k(i:end).*(Sv(i:end)-Su(i:end)))/k(i))/2;
  Kpsi(i) = (Sv(i)+trapz(log(k(i:end)), k(i:end).*(Sv(i:end)-Su(i:end)))/k(i))/2;
end

% Hydrostatic wave-vortex decomposition
% This decomposes the total energy E = 1/2(Su+Sv+Sb) into a wave component Ew
% and a residual vortical component Ev assuming hydrostatic waves. [uncomment
% for use]
%E = (Su+Sv+Sb)/2;
%Ew = 2*Kphi;
%Ev = E-Ew;

% Nonhydrostatic wave-vortex decomposition
% This decomposes the total energy E = 1/2(Su+Sv+Sw+Sb) into a wave component Ew
% and a residual vortical component Ev assuming hydrostatic waves.
E = (Su+Sv+Sw+Sb)/2;
Ew = 2*Kphi+Sw;
Ev = E-Ew;

% Plot the Helmholtz decomposition result.
figure()
hold on
loglog(k, K)
loglog(k, Kpsi)
loglog(k, Kphi)
xlabel('wavenumber [rad/m]')
ylabel('power spectral density [kg/s^2]')
legend('horizontal kinetic energy', 'rotational component', 'divergent component')

% Plot the wave-vortex decomposition result.
figure()
hold on
loglog(k, E)
loglog(k, Ew)
loglog(k, Ev)
xlabel('wavenumber [rad/m]')
ylabel('power spectral density [kg/s^2]')
legend('total energy', 'diagnosed wave component', 'residual vortical component')
