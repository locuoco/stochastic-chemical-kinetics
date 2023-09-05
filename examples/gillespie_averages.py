'''
    Stochastic enzyme kinetics: Stochastic enzyme kinetics averages plot
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
import seaborn as sns
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

dPs = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
ts  = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
n_simulations = {'Exact': 1000, 'tQSSA': 10000, 'sQSSA': 10000}
init_conditions =  {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}
max_t = 9

for s in sim:
	for _ in range(n_simulations[s]):
		g[s].x = init_conditions[s]
		g[s].t = 0
		states = g[s].simulate(t_final=max_t)
		prods = [state.x[g[s].species.P] for state in states]
		dprods = np.diff(prods)
		times = [state.t for state in states]
		dPs[s] += dprods.tolist()
		ts[s] += times[1:]

ndiv = max_t*50 + 1
bins = np.linspace(0, max_t, ndiv)

dy = {}
for s in sim:
	dprods = np.asarray(dPs[s]) / n_simulations[s]
	dy[s], _ = np.histogram(ts[s], bins, weights=dprods)

y = {}
var = {}
dx = bins[1] - bins[0]
for s in sim:
	y[s] = np.cumsum(dy[s])
	var[s] = y[s]*dy[s]/dx

x = (bins[1:] + bins[:-1])/2

for s in sim:
	error = np.sqrt(var[s])

	plt.plot(x, y[s], label=s)
	plt.fill_between(x, 0, error, alpha=.3)

plt.ylim(0)
plt.xlabel('Time')
plt.ylabel('Products count (average)')
plt.legend()

plt.show()























