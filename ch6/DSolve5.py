#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/21 00:44

@author: Niwatori
"""

import numpy as np
from math import pi, sin
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


'''
# Runge-Kutta Methods of Order 1 through 4
# Input: Function f(x, y), stride h, initial value (x0, y0), chosen order
# Output: Numerical solution of y(x) at discrete points
'''

def RungeKutta(f, h, x0, y0, order = 4):
    res = np.ndarray([Nt, Nx])
    res[0] = y0

    # 1-stage-1-order
    if order == 1:
        for k in range(1, Nt):
            x = x0 + k * h
            y = res[k - 1]
            K1 = f(x, y)
            res[k] = y + h * K1

    # 2-stage-2-order
    if order == 2:
        for k in range(1, Nt):
            x = x0 + k * h
            y = res[k - 1]
            K1 = f(x, y)
            K2 = f(x + h, y + h * K1)
            res[k] = y + (h / 2) * (K1 + K2)

    # 3-stage-3-order
    if order == 3:
        for k in range(1, Nt):
            x = x0 + k * h
            y = res[k - 1]
            K1 = f(x, y)
            K2 = f(x + h / 2, y + (h / 2) * K1)
            K3 = f(x + h, y - h * K1 + 2 * h * K2)
            res[k] = y + (h / 6) * (K1 + 4 * K2 + K3)

    # 4-stage-4-order
    if order == 4:
        for k in range(1, Nt):
            x = x0 + k * h
            y = res[k - 1]
            K1 = f(x, y)
            K2 = f(x + h / 2, y + (h / 2) * K1)
            K3 = f(x + h / 2, y + (h / 2) * K2)
            K4 = f(x + h, y + h * K3)
            res[k] = y + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)

    return res


''' Numerical solution for PDE '''
X = 2 * pi      # Max x interval [0, X]
Nx = 20         # Number of x subintervals
T = 30          # Max t interval [0, T]
Nt = 10 * T     # Number of t subintervals
hx = X / Nx     # Stride for x
ht = 0.1        # Stride for r

A = np.zeros([Nx, Nx])
y0 = np.zeros(Nx)
for i in range(Nx):
    A[i, (i + 1) % Nx] = -1
    A[(i + 1) % Nx, i] = 1
    y0[i] = sin(i * hx)
f = lambda x, y: np.dot(A, y) / (2 * hx)
# u = RungeKutta(f, ht, 0, y0, order = 1)
# u = RungeKutta(f, ht, 0, y0, order = 2)
# u = RungeKutta(f, ht, 0, y0, order = 3)
# u = RungeKutta(f, ht, 0, y0, order = 4)

''' Plot numerical results '''
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
t, x = np.mgrid[0:T-1e-8:ht, 0:X:hx]
ax.plot_surface(t, x, u, cmap = plt.cm.coolwarm, alpha = 0.4)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('u')
# plt.show()
# plt.savefig('fig5-4.svg', format = 'svg')