'''
    Stochastic enzyme kinetics: Stochastic enzyme kinetics averages plot using Gillespie
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

hists = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
n_simulations = 5000
init_conditions = {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}
max_t = 9

ndiv = max_t*50 + 1
bins = np.linspace(0, max_t, ndiv)

for s in sim:
	for _ in range(n_simulations):
		g[s].x = init_conditions[s]
		g[s].t = 0
		x, t = g[s].simulate(t_final=max_t)
		dprods = np.diff(x[:,g[s].species.P])
		dhist, _ = np.histogram(t[1:], bins, weights=dprods)
		hists[s].append(np.cumsum(dhist))
	hists[s] = np.array(hists[s])

avePs = {}
msqPs = {}
for s in sim:
	avePs[s] = np.sum(hists[s], axis=0) / hists[s].shape[0]
	msqPs[s] = np.sum(hists[s]**2, axis=0) / hists[s].shape[0]

t = (bins[1:] + bins[:-1])/2

for s in sim:
	error = np.sqrt(np.maximum(msqPs[s] - avePs[s]**2, 0))

	plt.plot(t, avePs[s], label=s)
	plt.fill_between(t, 0, error, alpha=.3)

plt.ylim(0)
plt.xlabel('Time')
plt.ylabel('Products count (average and st. dev., $P$)')
plt.legend()

plt.show()























