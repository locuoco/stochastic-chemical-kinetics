'''
    Stochastic enzyme kinetics: Stochastic enzyme kinetics averages plot using CME
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

avePs = {}
msqPs = {}
ts = {}
marginal_axes = {'Exact': (int(c['Exact'].species.C)+1), 'tQSSA': (), 'sQSSA': ()}
max_t = 9
Ps = np.arange(0, ST+1)

for s in sim:
	prob, ts[s] = c[s].simulate(dt=1e-4, t_final=max_t, n_sampling=100)
	marginal_prob = np.sum(prob, axis=marginal_axes[s])
	avePs[s] = np.tensordot(marginal_prob, Ps, axes=1)
	msqPs[s] = np.tensordot(marginal_prob, Ps**2, axes=1)

for s in sim:
	error = np.sqrt(np.maximum(msqPs[s] - avePs[s]**2, 0))

	plt.plot(ts[s], avePs[s], label=s)
	plt.fill_between(ts[s], 0, error, alpha=.3)

plt.ylim(0)
plt.xlabel('Time')
plt.ylabel('Products count (average and st. dev., $P$)')
plt.legend()

plt.show()























