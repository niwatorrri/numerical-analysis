#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/4/12 16:27

@author: Niwatori
"""

'''
# Solution to System of Nonlinear Equations (II)
#     x1 = -cos(x1) / 81 + x2 ** 2 / 9 + sin(x3) / 3
#     x2 = sin(x1) / 3 + cos(x3) / 3
#     x3 = -cos(x1) / 9 + x2 / 3 + sin(x3) / 6
# Output: Approximate solution
'''

import numpy as np
from math import sin, cos


''' Fixed-point iteration method '''

def SolveFPI(x):
    # steps = 0
    while True:
        # print('%d&%.15f\\\\'%(
        #     steps, np.linalg.norm(x - [0, 1 / 3, 0], 8)))
        # steps += 1
        x1 = np.array([-cos(x[0]) / 81 + x[1] ** 2 / 9 + sin(x[2]) / 3,
                       sin(x[0]) / 3 + cos(x[2]) / 3,
                       -cos(x[0]) / 9 + x[1] / 3 + sin(x[2]) / 6])
        if np.linalg.norm(x1 - x, 8) < 1e-15: break
        x = x1
    return x1


''' Newton's method '''

def SolveNewton(x):
    # steps = 0
    while True:
        # print('%d&%.15f\\\\'%(
        #     steps, np.linalg.norm(x - [0, 1 / 3, 0])))
        # steps += 1
        A = np.matrix([[sin(x[0]) / 81 - 1, x[1] * 2 / 9, cos(x[2]) / 3],
                       [cos(x[0]) / 3, -1, -sin(x[2]) / 3],
                       [sin(x[0]) / 9, 1 / 3, cos(x[2]) / 6 - 1]])
        b = np.array([-cos(x[0]) / 81 + x[1] ** 2 / 9 + sin(x[2]) / 3 - x[0],
                      sin(x[0]) / 3 + cos(x[2]) / 3 - x[1],
                      -cos(x[0]) / 9 + x[1] / 3 + sin(x[2]) / 6 - x[2]])
        y = np.linalg.solve(A, -b)
        if np.linalg.norm(y) < 1e-15: break
        x = x + y
    return x


print(SolveFPI(np.array([0, 0, 0])))
print(SolveNewton(np.array([0, 0, 0])))
