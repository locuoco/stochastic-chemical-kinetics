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

x, t = g.simulate()

plt.plot(t, x[:,g.species.C], drawstyle='steps-post', label='C')
plt.plot(t, x[:,g.species.P], drawstyle='steps-post', label='P')

plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()

plt.show()























