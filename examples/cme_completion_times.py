'''
    Stochastic enzyme kinetics: completion times histogram using CME
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
import sys

sys.path.append('../pybind')

import cme

kf = 10
kb = 9
kcat = 1
kM = (kb + kcat) / kf
ET = 10
ST = 9

sim = ['Exact', 'tQSSA', 'sQSSA']
c = {}

c['Exact'] = cme.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
c['tQSSA'] = cme.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
c['sQSSA'] = cme.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

bins = np.linspace(0, 9, 19)
completion_time_weights = {
	'Exact': np.empty_like(bins),
	'tQSSA': np.empty_like(bins),
	'sQSSA': np.empty_like(bins)
}
completion_states = {'Exact': (0, ST), 'tQSSA': ST, 'sQSSA': ST}

for s in sim:
	for bin, i in zip(bins, range(bins.size)):
		c[s].simulate(dt=1e-4, t_final=bin, noreturn=True)
		completion_time_weights[s][i] = c[s].p[completion_states[s]]

for s in sim:
	completion_time_weights[s] = np.diff(completion_time_weights[s])

def weighted_ave_sd(arr, weights):
	average = np.average(arr, weights=weights)
	variance = np.average((arr-average)**2, weights=weights)
	return average, np.sqrt(variance)

for s in sim:
	average, sd = weighted_ave_sd(bins[1:], weights=completion_time_weights[s])
	print(s, average, '+/-', sd)

t = (bins[1:] + bins[:-1])/2

plt.hist(t, bins, weights=completion_time_weights['Exact'], label='Exact', density=True, color='blue', histtype='step', fill=False)
plt.hist(t, bins, weights=completion_time_weights['tQSSA'], label='tQSSA', density=True, color='red', alpha=.6)
plt.hist(t, bins, weights=completion_time_weights['sQSSA'], label='sQSSA', density=True, color='gray', alpha=.3)

plt.ylim(0, .4)
plt.xlabel('Completion time (τ)')
plt.ylabel('Probability density')
plt.legend()

plt.show()























