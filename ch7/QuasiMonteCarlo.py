#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/6/8 00:08

@author: Niwatori
"""

import random
from math import pi, sin, log

'''
# Quasi-Monte-Carlo Method: Numerical integration using Halton sequence
# Input: Function f, Dimension d, Number of samples N
# Output: Approximate value of integral of f over [0,1]^d
'''

prime = [-1, 2, 3, 5, 7, 11, 13, 17, 19]

def Halton(n, k):
    p = prime[k]
    res = r = 0
    while n > 0:
        res = res * p + n % p
        n = n // p
        r = r + 1
    return res / (p ** r)


def QuasiMC(f, d, N):
    sum = 0
    for n in range(1, N + 1):
        x = []
        for k in range(1, d + 1):
            x.append(Halton(n, k))
        sum += f(x)
    return sum / N


def MonteCarlo(f, N):
    res = 0
    for k in range(100):
        I = 0
        for i in range(N):
            x = [random.random(), random.random()]
            I = I + f(x)
        res += I / N
    res = res / 100
    return res


f = lambda x: 4 * x[1] * sin(pi * x[0])
ans = 4 / pi
N = 1
err1 = []
err2 = []

for k in range(6):
    N = N * 10
    err1.append(abs(QuasiMC(f, 2, N) - ans))
    err2.append(abs(MonteCarlo(f, N) - ans))
    print('%d & %.6f & '%(N, err1[-1]), end = '')
    if k > 0:
        print('%.6f'%((log(err1[-2]) - log(err1[-1])) / log(10)), end = '')
    print(' & %.6f & '%(err2[-1]), end = '')
    if k > 0:
        print('%.6f'%((log(err2[-2]) - log(err2[-1])) / log(10)), end = '')
    print('\\\\')
