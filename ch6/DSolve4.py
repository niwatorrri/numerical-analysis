#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/24 19:58

@author: Niwatori
"""

import numpy as np
from math import log, sin, cos, copysign
sgn = lambda x: copysign(1, x)

'''
# Shooting Method for Boundary Value Problems
# Input: Function f of equation system, initial guess s in [a, b] and tolerance
# Output: Numerical solution for BVP
'''

def Shoot(f, a, b):
    h = 1e-2

    ''' Bisection method to give a feasible initial guess of s '''
    # def phi(s):
    #     y = RungeKutta(f, h, 1, np.array([1, s, 0, 1]))[-1]
    #     return y[0] - 2
    #
    # fa = phi(a)
    # fb = phi(b)
    # if sgn(fa) * sgn(fb) > 0:
    #     return 'f(a)f(b)<0 not satisfied.'
    #
    # while (b - a) / 2 > 1e-1:
    #     c = (a + b) / 2
    #     fc = phi(c)
    #     if abs(fc) < 1e-8:
    #         break
    #     if sgn(fa) * sgn(fc) < 0:
    #         b = c; fb = fc
    #     else:
    #         a = c; fa = fc
    s = (a + b) / 2

    ''' Newton's method for faster convergence '''
    while True:
        y = RungeKutta(f, h, 1, np.array([1, s, 0, 1]))[-1]
        delta = (y[0] - 2) / y[2]
        s = s - delta
        if abs(delta) < 1e-8:
            break
    res = RungeKutta(f, h, 1, np.array([1, s, 0, 1]))
    return [t[0] for t in res]


''' Classical Runge-Kutta Method '''

def RungeKutta(f, h, x0, y0):
    x = x0
    y = y0
    res = [y0]

    while x + 1e-8 < 2:
        K1 = f(x, y)
        K2 = f(x + h / 2, y + (h / 2) * K1)
        K3 = f(x + h / 2, y + (h / 2) * K2)
        K4 = f(x + h, y + h * K3)
        y = y + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        x = x + h
        res.append(y)
    return res


''' System of differential equations '''
def f(x, y):
    f1 = y[1]
    f2 = (sin(log(x)) + 2 * y[0]) / (x ** 2) - (2 * y[1] / x)
    f3 = y[3]
    f4 = (2 * y[2]) / (x ** 2) - (2 * y[3] / x)
    return np.array([f1, f2, f3, f4])


''' Accurate solution '''
c2 = (8 - 12 * sin(log(2)) - 4 * cos(log(2))) / 70
c1 = (11 / 10) - c2
g = lambda x: (c1 * x) + c2 / (x ** 2) - 3 * sin(log(x)) / 10 - cos(log(x)) / 10

''' Numerical solution '''
y = Shoot(f, -5, 5)

for i in range(1, 11):
    x = 1 + i * 0.1
    num = y[i * 10]
    acc = g(x)
    print('%.1f&%.12f&%.12f&%.2e\\\\'%(x, num, acc, abs(num - acc)))
