#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/11 19:15

@author: Niwatori
"""

from math import pi, tan, log

'''
# Euler's Method (Forward / Predictor-Corrector)
# Input: Function f(x, y), stride h, initial value (x0, y0), chosen mode
# Output: Numerical solution of y(x) at discrete points
'''

def Euler(f, h, x0, y0, mode = "forward"):
    x = x0
    y = y0
    res = [[x0, y0]]

    # Forward Euler's method
    if mode == "forward":
        while x < 2:
            y = y + h * f(x, y)
            x = x + h
            res.append([x, y])

    # Modified Euler's method - Predictor-Corrector method
    elif mode == "modified":
        while x < 2:
            K1 = f(x, y)
            K2 = f(x + h, y + h * K1)
            y = y + (K1 + K2) * h / 2
            x = x + h
            res.append([x, y])

    return res


f = lambda x, y: -1 / (x ** 2) - (y / x) - (y ** 2)
g = lambda x: -tan(log(x) + pi / 4) / x     # Exact solution
h = 0.01

# Numerical solution
res1 = Euler(f, h, x0 = 1, y0 = -1, mode = "forward")
res2 = Euler(f, h, x0 = 1, y0 = -1, mode = "modified")
res3 = [[t / 100, g(t / 100)] for t in range(100, 201)]
for i in range(1, 11):
    print("%.1f&%.6f&%.6f&%.6f\\\\"%(
        1 + i / 10, res1[i * 10][1], res2[i * 10][1], res3[i * 10][1]))
