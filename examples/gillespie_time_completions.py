'''
    Stochastic enzyme kinetics: Stochastic enzyme kinetics time completions histogram
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

g_ss = gillespie.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
g_tq = gillespie.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
g_sq = gillespie.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

completion_times_ss = []
completion_times_tq = []
completion_times_sq = []

n_simulations = 20000
max_t = 20

for _ in range(n_simulations):
	g_ss.x = [0, 0]
	g_ss.t = 0
	g_ss.simulate(t_final=max_t, noreturn=True)
	completion_times_ss.append(g_ss.t)

for _ in range(n_simulations):
	g_tq.x = [0]
	g_tq.t = 0
	g_tq.simulate(t_final=max_t, noreturn=True)
	completion_times_tq.append(g_tq.t)

for _ in range(n_simulations):
	g_sq.x = [0]
	g_sq.t = 0
	g_sq.simulate(t_final=max_t, noreturn=True)
	completion_times_sq.append(g_sq.t)

print('Exact:', mean(completion_times_ss), '+/-', stdev(completion_times_ss))
print('tQSSA:', mean(completion_times_tq), '+/-', stdev(completion_times_tq))
print('sQSSA:', mean(completion_times_sq), '+/-', stdev(completion_times_sq))

bins = np.linspace(0, 9, 19)

plt.ylim(0, .4)
plt.hist(completion_times_ss, bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
plt.hist(completion_times_tq, bins, label='tQSSA', density=True, color='red', alpha=.6)
plt.hist(completion_times_sq, bins, label='sQSSA', density=True, color='gray', alpha=.3)

plt.xlabel('Time')
plt.ylabel('Probability density')
plt.legend()

plt.show()























