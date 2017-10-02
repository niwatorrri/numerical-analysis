#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/4/12 20:21

@author: Niwatori
"""

'''
# Solution to System of Nonlinear Equations (III)
# Both use Newton's method
# Output: Approximate solution
'''

import numpy as np
from math import sqrt, exp


'''
# Equation 3:
#    x1^2 + x2^2 + x3^2 = 5
#    x1 + x2 = 1
#    x1 + x3 = 3
'''

def Solve1(x):
    # steps = 0
    while True:
        # print('%d&(%.4e, %.4e, %.4e)&%.4e\\\\' % (
        #     steps, x[0], x[1], x[2], np.linalg.norm(x - [5 / 3, -2 / 3, 4 / 3])))
        # steps += 1
        A = np.matrix([[2 * x[0], 2 * x[1], 2 * x[2]], [1, 1, 0], [1, 0, 1]])
        b = np.array([x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 5, x[0] + x[1] - 1, x[0] + x[2] - 3])
        y = np.linalg.solve(A, -b)
        if np.linalg.norm(y) < 1e-15: break
        x = x + y
    return x


'''
# Equation 4:
#    10000 * x1 * x2 = 1
#    exp(-x1) + exp(-x2) = 1.0001
'''

def Solve2(x):
    # steps = 0
    while True:
        # print('%d&(%.4e, %.4e)&%.4e\\\\' % (
        #     steps, x[0], x[1], np.linalg.norm(x - [ 1.09815933e-05, 9.10614674e+00])))
        # steps += 1
        A = np.matrix([[10000 * x[1], 10000 * x[0]], [-exp(-x[0]), -exp(-x[1])]])
        b = np.array([10000 * x[0] * x[1] - 1, exp(-x[0]) + exp(-x[1]) - 1.0001])
        y = np.linalg.solve(A, -b)
        if np.linalg.norm(y) < 1e-15: break
        x = x + y
    return x


print(Solve1(np.array([(1 + sqrt(3)) / 2, (1 - sqrt(3)) / 2, sqrt(3)])))
print(Solve2(np.array([0, 1])))
# print(Solve2(np.array([1, 0])))
