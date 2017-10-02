#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/19 21:34

@author: Niwatori
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


'''
# Classical Runge-Kutta Method
# Input: Function f(x, y), stride h, initial value (x0, y0)
# Output: Numerical solution of y(x) at discrete points
'''

def RungeKutta(f, h, x0, y0):
    x = x0
    y = y0
    res = [y0]

    while x + 1e-10 < 50:
        K1 = f(x, y)
        K2 = f(x + h / 2, y + (h / 2) * K1)
        K3 = f(x + h / 2, y + (h / 2) * K2)
        K4 = f(x + h, y + h * K3)
        y = y + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        x = x + h
        res.append(y)
    return res


''' Lorenz System '''
def func(sigma, rho, beta):
    def f(t, u):
        f1 = sigma * (u[1] - u[0])
        f2 = u[0] * (rho - u[2]) - u[1]
        f3 = u[0] * u[1] - beta * u[2]
        return np.array([f1, f2, f3])
    return f


''' Numerical solution with different initial values '''
f = func(sigma = 10, rho = 28, beta = 8 / 3)
h = 0.001
# res1 = RungeKutta(f, h, 0, np.array([1, 1, 1]))
# res2 = RungeKutta(f, h, 0, np.array([1 + 1e-5, 1, 1]))
# res3 = RungeKutta(f, h, 0, np.array([100, 100, 100]))

g = func(sigma = 10, rho = 0.5, beta = 8 / 3)
# res4 = RungeKutta(g, h, 0, np.array([5, 5, 5]))
# res5 = RungeKutta(g, h, 0, np.array([10, 10, 10]))
# res6 = RungeKutta(g, h, 0, np.array([20, 20, 20]))

g = func(sigma = 10, rho = 14, beta = 8 / 3)
# res4 = RungeKutta(g, h, 0, np.array([0, 0, 20]))
# res5 = RungeKutta(g, h, 0, np.array([20, 0, 0]))
# res6 = RungeKutta(g, h, 0, np.array([15, 15, 15]))
# res7 = RungeKutta(g, h, 0, np.array([20, 20, 20]))


''' 3D Parametric plot '''
fig = plt.figure()
ax = fig.gca(projection = '3d')
x = lambda res: [t[0] for t in res]
y = lambda res: [t[1] for t in res]
z = lambda res: [t[2] for t in res]
# ax.plot(x(res1), y(res1), z(res1), label = '(1, 1, 1)')
# ax.plot(x(res2), y(res2), z(res2), label = '(1+1e-5, 1, 1)')
# ax.plot(x(res3), y(res3), z(res3), label = '(100, 100, 100)')
# ax.plot(x(res4), y(res4), z(res4), label = '(0, 0, 20)')
# ax.plot(x(res5), y(res5), z(res5), label = '(20, 0, 0)')
# ax.plot(x(res6), y(res6), z(res6), label = '(15, 15, 15)')
# ax.plot(x(res7), y(res7), z(res7), label = '(20, 20, 20)')
ax.legend()
# plt.show()
# plt.savefig("fig8.eps")