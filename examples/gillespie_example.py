'''
    Stochastic enzyme kinetics: Gillespie algorithm python example
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

g = gillespie.single_substrate(kf=10, kb=9, kcat=1, ET=10, ST=9)

states = g.simulate()

Cs = [state.x[g.species.C] for state in states]
Ps = [state.x[g.species.P] for state in states]
times = [state.t for state in states]

plt.plot(times, Cs, drawstyle='steps-post', label='C')
plt.plot(times, Ps, drawstyle='steps-post', label='P')

plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()

plt.show()























