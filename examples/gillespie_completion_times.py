'''
    Stochastic enzyme kinetics: completion times histogram
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
import sys

sys.path.append('../pybind')

import gillespie

kf = 10
kb = 9
kcat = 1
kM = (kb + kcat) / kf
ET = 10
ST = 9

sim = ['Exact', 'tQSSA', 'sQSSA']
g = {}

g['Exact'] = gillespie.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
g['tQSSA'] = gillespie.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
g['sQSSA'] = gillespie.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

completion_times = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
init_conditions =  {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}

n_simulations = 20000
max_t = 20

for s in sim:
	for _ in range(n_simulations):
		g[s].x = init_conditions[s]
		g[s].t = 0
		g[s].simulate(t_final=max_t, noreturn=True)
		completion_times[s].append(g[s].t)

for s in sim:
	print(s, mean(completion_times[s]), '+/-', stdev(completion_times[s]))

bins = np.linspace(0, 9, 19)

plt.hist(completion_times['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
plt.hist(completion_times['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
plt.hist(completion_times['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

plt.ylim(0, .4)
plt.xlabel('Completion time ($\tau$)')
plt.ylabel('Probability density')
plt.legend()

plt.show()























