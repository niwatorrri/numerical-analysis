#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/4/9 23:39

@author: Niwatori
"""

'''
# Solution to System of Nonlinear Equations (I)
#     (x_1 + 3) * (x_2 ^ 3 - 7) + 18 = 0
#     sin(x_2 * exp(x_1) - 1) = 0
# Output: Approximate solution
'''

import numpy as np
from math import exp, sin, cos


''' Newton's method '''
''' A is the Jacobian matrix '''

def SolveNewton(x):
    # steps = 0
    while True:
        # print('%d&(%.10f, %.10f)&%.10f\\\\'%(
        #     steps, x[0], x[1], np.linalg.norm(x - [0, 1])))
        # steps += 1
        A = np.matrix([[x[1] ** 3 - 7, 3 * (x[0] + 3) * x[1] ** 2],
                       [x[1] * exp(x[0]) * cos(x[1] * exp(x[0]) - 1), exp(x[0]) * cos(x[1] * exp(x[0]) - 1)]])
        b = np.array([(x[0] + 3) * (x[1] ** 3 - 7) + 18, sin(x[1] * exp(x[0]) - 1)])
        y = np.linalg.solve(A, -b)
        if np.linalg.norm(y) < 1e-15: break
        x = x + y
    return x


''' Broyden's method '''
''' A is the approximate inverse of Jacobian matrix'''

def SolveBroyden(x):
    f = lambda x: \
        np.matrix([(x[0, 0] + 3) * (x[1, 0] ** 3 - 7) + 18, sin(x[1, 0] * exp(x[0, 0]) - 1)]).T
    A = np.matrix([[x[1] ** 3 - 7, 3 * (x[0] + 3) * x[1] ** 2],
                   [x[1] * exp(x[0]) * cos(x[1] * exp(x[0]) - 1), exp(x[0]) * cos(x[1] * exp(x[0]) - 1)]]).I
    x = np.matrix(x).T
    x1 = x - A * f(x)
    # steps = 0

    while np.linalg.norm(x1 - x) > 1e-15:
        # print('%d&(%.10f, %.10f)&%.10f\\\\'%(
        #     steps, x[0], x[1], np.linalg.norm(x - np.matrix([0, 1]).T)))
        # steps += 1
        y = x1 - x
        g = f(x1) - f(x)
        A = A - ((A * g - y) * y.T * A) / (y.T * A * g)
        x = x1
        x1 = x1 - A * f(x1)
    return x1


print(SolveNewton(np.array([-0.5, 1.4])))
print(SolveBroyden(np.array([-0.5, 1.4])))
