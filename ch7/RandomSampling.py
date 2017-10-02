#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/25 23:35

@author: Niwatori
"""

import random
from math import pi, sin, cos, acos, log
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
# Generate a random vector with p.d.f.
#   p(x1, x2, x3) = exp(-r) / (8pi)
'''

x = []
y = []
z = []
N = 200
for t in range(N):
    u1 = random.random()
    u2 = random.random()
    u3 = random.random()
    u4 = random.random()
    u5 = random.random()
    r = -log(u1) - log(u2) - log(u3)
    phi = acos(1 - 2 * u4)
    theta = 2 * pi * u5
    x.append(r * sin(phi) * cos(theta))
    y.append(r * sin(phi) * sin(theta))
    z.append(r * cos(phi))

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(x, y, z)
plt.savefig('fig.eps')
