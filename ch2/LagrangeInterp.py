#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 23:30:52 2017

@author: Niwatori
"""

'''
# Lagrange's Interpolation method
# Input: Data points, corresponding values and number of points
# Output: Interpolation polynomial
'''
import pylab as pl
from math import pi, cos

def LagrangeInterp(x, y, n):
    def InterpPoly(t):
        ans = 0
        for i in range(n):
            tmp = 1
            for j in range(n):
                if j != i:
                    tmp = tmp * (t-x[j]) / (x[i]-x[j])
            ans = ans + y[i] * tmp
        return ans
    return InterpPoly

# Runge function r(x) and approximate polynomial p(x)
x = [5*cos(t*pi/42) for t in range(1, 42, 2)]
r = lambda x:1/(1+x**2)
p = LagrangeInterp(x, [r(t) for t in x], 21)

# Plot images of two functions
X = pl.linspace(-5, 5, 100)
pl.plot(X, [r(t) for t in X], label="Runge", linewidth=2)
pl.plot(X, [p(t) for t in X], label="Approximation 2")
pl.legend(loc='upper right')
pl.savefig("fig2.eps")