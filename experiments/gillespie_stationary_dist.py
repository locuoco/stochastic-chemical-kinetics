'''
    Stochastic enzyme kinetics: Goldbeter-Koshland switch stationary
        distribution using Gillespie algorithm
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

import gillespie

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

sim = ['Exact', 'tQSSA', 'sQSSA']
g = {}

g['Exact'] = gillespie.goldbeter_koshland(kfe=kfe, kbe=kbe, ke=ke, kfd=kfd, kbd=kbd, kd=kd, ET=ET, DT=DT, ST=ST)
g['tQSSA'] = gillespie.goldbeter_koshland_tqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)
g['sQSSA'] = gillespie.goldbeter_koshland_sqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)

SP_hats = {}
init_conditions = {'Exact': [ST//2, 0, 0], 'tQSSA': [ST//2], 'sQSSA': [ST//2]}
max_t = 1000

tracked_species = {
	'Exact': np.array([g['Exact'].species.SP, g['Exact'].species.CP], dtype=int),
	'tQSSA': np.array([g['tQSSA'].species.SP_hat], dtype=int),
	'sQSSA': np.array([g['sQSSA'].species.SP], dtype=int)
}

for s in sim:
	g[s].x = init_conditions[s]
	g[s].t = 0
	x, _ = g[s].simulate(t_final=max_t)
	SP_hats[s] = np.sum(x[:,tracked_species[s]], axis=1)

for s in sim:
	print(s, np.mean(SP_hats[s]), '+/-', np.std(SP_hats[s]))

bins = np.linspace(0, ST, 21)

plt.hist(SP_hats['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
plt.hist(SP_hats['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
plt.hist(SP_hats['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

plt.xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
plt.ylabel('Probability density')
plt.legend()

plt.show()























