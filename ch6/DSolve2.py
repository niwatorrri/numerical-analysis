#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/18 18:14

@author: Niwatori
"""

import numpy as np
import pylab as pl
import time


'''
# Runge-Kutta Methods: Modified Euler (2) and Classical Runge-Kutta (4)
# Input: Function f(x, y), stride h, initial value (x0, y0), chosen order
# Output: Numerical solution of y(x) at discrete points
'''

def RungeKutta(f, h, x0, y0, order = 4):
    x = x0
    y = y0
    res = [y0]

    # 2-stage-2-order (Modified Euler)
    if order == 2:
        while x + 1e-10 < 10:
            K1 = f(x, y)
            K2 = f(x + h, y + h * K1)
            y = y + (h / 2) * (K1 + K2)
            x = x + h
            res.append(y)

    # 4-stage-4-order (Classical Runge-Kutta)
    if order == 4:
        while x + 1e-10 < 10:
            K1 = f(x, y)
            K2 = f(x + h / 2, y + (h / 2) * K1)
            K3 = f(x + h / 2, y + (h / 2) * K2)
            K4 = f(x + h, y + h * K3)
            y = y + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
            x = x + h
            res.append(y)

    return res


'''
# Predictor-Corrector Method: Adams-Bashforth Method (Four-point)
# Input: Function f(x, y), stride h, initial value (x0, y0)
# Output: Numerical solution of y(x) at discrete points
'''

def Adams4(f, h, x0, y0):
    x = x0
    y = [y0]
    n = 0

    # Initial four points (with modified Euler)
    while n < 3:
        K1 = f(x, y[n])
        K2 = f(x + h, y[n] + h * K1)
        ynew = y[n] + (h / 2) * (K1 + K2)
        x = x + h
        n = n + 1
        y.append(ynew)

    # Four-point Adams-Bashforth predictor-corrector method
    while x + 1e-10 < 10:
        u = y[n] + h * (55 * f(x, y[n]) - 59 * f(x - h, y[n - 1])
                        + 37 * f(x - 2 * h, y[n - 2]) - 9 * f(x - 3 * h, y[n - 3])) / 24
        ynew = y[n] + h * (9 * f(x + h, u) + 19 * f(x, y[n])
                        - 5 * f(x - h, y[n - 1]) + f(x - 2 * h, y[n - 2])) / 24
        x = x + h
        n = n + 1
        y.append(ynew)

    return y


''' The function to be solved '''
def f(x, y):
    f1 = -0.013 * y[0] - 1000 * y[0] * y[1]
    f2 = -2500 * y[1] * y[2]
    f3 = -0.013 * y[0] - 1000 * y[0] * y[1] - 2500 * y[1] * y[2]
    return np.array([f1, f2, f3])


''' Numerical solutions '''
h = 0.0001
# start = time.clock()
res1 = RungeKutta(f, h, 0, np.array([1, 1, 0]), order = 2)
res2 = RungeKutta(f, h, 0, np.array([1, 1, 0]), order = 4)
res3 = Adams4(f, h, 0, np.array([1, 1, 0]))
# end = time.clock()
# print(end - start)


''' Print solution within [0, 0.005] '''
for i in range(1, 11):
    for j in range(0, 3):
        if j == 1:
            print('%.4f'%(i * 0.0005), end = '')
        print('&%.6f&%.6f&%.6f\\\\'%(res1[i * 5][j], res2[i * 5][j], res3[i * 5][j]))
    print('\hline')


''' Plot the integral curves '''
X = pl.linspace(0, 0.01, num = 101)
pl.plot(X, [res3[t][0] for t in range(0, 101)], label = "y1(t)")
pl.plot(X, [res3[t][1] for t in range(0, 101)], label = "y2(t)")
pl.plot(X, [res3[t][2] for t in range(0, 101)], label = "y3(t)")
pl.legend(loc = 'center right')
pl.xlim(0, 0.01)
pl.ylim(-0.5, 2.3)
pl.savefig("fig1.eps")