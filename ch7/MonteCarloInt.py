#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/25 22:04

@author: Niwatori
"""

import random
from math import cos, exp, log

'''
# Monte Carlo Methods for Numerical Integration
# View an integral as the expectation of a random variable
'''


''' Problem A '''
N = 1
res = []
ans = 1 - exp(-2)

for t in range(5):              # For each sample size
    N = N * 10
    aveI = err = 0
    for k in range(100):        # Repeat 100 times
        I = 0
        for i in range(N):      # Sampling
            x = random.random()
            I = I + 2 * exp(-2 * x)
        aveI = aveI + I / N
        err = err + abs(I / N - ans)
    aveI = aveI / 100           # Average integral
    err = err / 100             # Average error
    res.append([N, aveI, err])

for i in range(5):
    print('%d&%.6f&%.6f&'%(res[i][0], res[i][1], res[i][2]), end = '')
    if i > 0:
        print('%.4f'%((log(res[i - 1][2]) - log(res[i][2])) / log(10)), end = '')
    print('\\\\')


''' Problem B '''
N = 1
res1 = []
res2 = []
ans = 0.198373
C = 5 / (1 - exp(-5))

for t in range(4):              # For each sample size
    N = N * 10
    ave1 = ave2 = err1 = err2 = 0
    for k in range(100):        # Repeat 100 times
        I1 = I2 = 0
        for i in range(N):      # Sampling
            x = random.random()
            I1 = I1 + cos(x / 5) * exp(-5 * x)
            y = -log(1 - 5 * x / C) / 5
            I2 = I2 + cos(y / 5) / C
        ave1 += I1 / N
        ave2 += I2 / N
        err1 += abs(I1 / N - ans)
        err2 += abs(I2 / N - ans)
    res1.append([N, ave1 / 100, err1 / 100])    # Average integral
    res2.append([N, ave2 / 100, err2 / 100])    # Average error

for i in range(4):
    print('%d&%.6f&%.6f&%.6f&%.6f\\\\'%(res1[i][0],
            res1[i][1], res2[i][1], res1[i][2], res2[i][2]))
