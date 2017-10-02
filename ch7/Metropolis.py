#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2017/5/30 13:42

@author: Niwatori
"""

import numpy, random
import pylab as pl
from math import exp

'''
# Metropolis-Hastings Algorithm (2-D Ising Model)
# Markov Chain Monte Carlo Method
'''

def Solve(M, beta):
    N = 10000
    s = numpy.ones((M, M))
    h = H(s)
    U = C = 0
    for k in range(N // 5 * 4):
        x, y, ok, delta = Metropolis(s, beta)
        if ok == 1: s[x][y] = -s[x][y]
        h = h + delta
    for k in range(N // 5):
        x, y, ok, delta = Metropolis(s, beta)
        if ok == 1: s[x][y] = -s[x][y]
        h = h + delta
        U = U + h
        C = C + h ** 2
    U = U / N * 5
    C = C / N * 5
    return [U / (M ** 2), (C - U ** 2) * (beta / M) ** 2]


def Metropolis(s, beta):
    x = random.randint(0, M - 1)
    y = random.randint(0, M - 1)
    delta = 0
    delta += 2 * (s[x][y] * s[(x - 1 + M) % M][y])
    delta += 2 * (s[x][y] * s[(x + 1 + M) % M][y])
    delta += 2 * (s[x][y] * s[x][(y - 1 + M) % M])
    delta += 2 * (s[x][y] * s[x][(y + 1 + M) % M])
    if random.random() < min(1, exp(-beta * delta)):
        return [x, y, 1, delta]
    else:
        return [x, y, 0, 0]


def H(s):
    sum = 0
    for i in range(M):
        for j in range(M):
            sum += s[i][j] * s[(i + 1) % M][j]
            sum += s[i][j] * s[i][(j + 1) % M]
    return -sum


M = 200
x = pl.linspace(0, 1, num = 101)
res = [Solve(M, beta) for beta in x]
pl.plot(x, [t[0] for t in res], label = "U_M")
pl.plot(x, [t[1] for t in res], label = "C_M", color = "green")

pl.legend(loc = 'upper right')
pl.show()
# pl.savefig("figU.eps")