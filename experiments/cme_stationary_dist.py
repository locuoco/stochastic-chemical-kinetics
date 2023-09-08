'''
    Stochastic enzyme kinetics: Goldbeter-Koshland switch stationary
        distribution using CME
    Copyright (C) 2023 Alessandro Lo Cuoco (alessandro.locuoco@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
import matplotlib.pyplot as plt
from math import ceil
import sys

sys.path.append('../pybind')

import cme

kfe = 10
kbe = 8.3
ke = 1.7
kfd = 10
kbd = 8.3
kd = 1.7
kME = (kbe + ke) / kfe
kMD = (kbd + kd) / kfd
ET = 100
DT = 100
ST = 100

# not using exact formulation due to very long computational time
sim = ['tQSSA', 'sQSSA']
c = {}

c['tQSSA'] = cme.goldbeter_koshland_tqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)
c['sQSSA'] = cme.goldbeter_koshland_sqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)

SP_hat_dists = {}

init_conditions = {}
for s in sim:
	init_conditions[s] = np.zeros_like(c[s].p)
init_conditions['tQSSA'][ST//2] = 1
init_conditions['sQSSA'][ST//2] = 1

max_t = 10
ave = {}
msq = {}
possible_SP_hats = np.arange(0, ST+1)

for s in sim:
	c[s].p = init_conditions[s]
	c[s].t = 0
	c[s].simulate(dt=1e-3, t_final=max_t, noreturn=True)
	SP_hat_dists[s] = c[s].p
	ave[s] = np.dot(c[s].p, possible_SP_hats)
	msq[s] = np.dot(c[s].p, possible_SP_hats**2)

for s in sim:
	print(s, ave[s], '+/-', np.sqrt(np.maximum(msq[s] - ave[s]**2, 0)))

bins = np.linspace(0, ST, 21)
bins_all = np.linspace(0, ST+1, ST+2)
x = (bins_all[1:] + bins_all[:-1])/2

plt.hist(x, bins, weights=SP_hat_dists['tQSSA'], label='tQSSA', density=True, color='red', alpha=.6)
plt.hist(x, bins, weights=SP_hat_dists['sQSSA'], label='sQSSA', density=True, color='gray', alpha=.3)

plt.xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
plt.ylabel('Probability density')
plt.legend()

plt.show()























