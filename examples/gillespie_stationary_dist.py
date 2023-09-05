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
from statistics import mean, stdev
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

SP_hats = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
init_conditions = {'Exact': [ST//2, 0, 0], 'tQSSA': [ST//2], 'sQSSA': [ST//2]}
n_samplings = {'Exact': ceil((kfe+kbe+kfd+kbd)/(ke+kd)), 'tQSSA': 1, 'sQSSA': 1}
max_t = 500

tracked_species = {
	'Exact': [g['Exact'].species.SP, g['Exact'].species.CP],
	'tQSSA': [g['tQSSA'].species.SP_hat],
	'sQSSA': [g['sQSSA'].species.SP]
}

for s in sim:
	g[s].x = init_conditions[s]
	g[s].t = 0
	states = g[s].simulate(t_final=max_t, n_sampling=n_samplings[s])
	SP_hats[s] = [sum([state.x[species] for species in tracked_species[s]]) for state in states]

for s in sim:
	print(s, mean(SP_hats[s]), '+/-', stdev(SP_hats[s]))

bins = np.linspace(0, ST, 21)

plt.hist(SP_hats['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
plt.hist(SP_hats['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
plt.hist(SP_hats['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

plt.xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
plt.ylabel('Probability density')
plt.legend()

plt.show()























