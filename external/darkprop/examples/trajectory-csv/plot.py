#!/usr/bin/env python3
"""
plot trajectories
"""

import os
import numpy as np
import matplotlib.pyplot as plt

BASEDIR = 'out'

if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for file in os.listdir(BASEDIR):
        x, y, z = np.loadtxt(os.path.join(BASEDIR, file), unpack=True)
        z -= 6371  # Earth's radius
        ax.scatter(x, y, z)
        ax.plot(x, y, z, ls='-')
        ax.set_xlabel('x [km]')
        ax.set_ylabel('y [km]')
        ax.set_zlabel('z [km]')

plt.show()
