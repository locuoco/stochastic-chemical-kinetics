'''
    Stochastic enzyme kinetics: CME python example
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

# have mmpeg installed to run this script. Using conda:

#   conda install ffmpeg

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, colors
import sys

sys.path.append('../pybind')

import cme

c = cme.single_substrate(kf=10, kb=9, kcat=1, ET=10, ST=9)
init_p = c.p.copy()

fig = plt.figure()
ax = plt.axes()
im = ax.imshow(c.p, norm=colors.SymLogNorm(1e-20), origin='lower')
plt.colorbar(im, ax=ax, label='Probability')
im.set_clim(1e-10, 1)

annot = ax.annotate('Time: 0.000', (0.09, 0.92), xycoords='figure fraction')

def init():
	c.p = init_p
	c.t = 0
	im.set_data(c.p)
	annot.set_text('Time: 0.000')
	return im, annot

def animate(i):
	current_t = i*0.002
	c.simulate(dt=1e-4, t_final=current_t, noreturn=True)
	im.set_data(c.p)
	annot.set_text('Time: {:.3f}'.format(current_t))
	return im, annot

fps = 30

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=5000, interval=1000/fps, blit=False)

if len(sys.argv) > 1:
	if sys.argv[1] == 's' or sys.argv[1] == 'save':
		anim.save('resources/cme_example.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

plt.xlabel('Product count ($P$)')
plt.ylabel('Enzyme-substrate complex count ($C$)')
plt.show()























