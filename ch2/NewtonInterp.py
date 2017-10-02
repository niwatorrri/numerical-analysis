#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 20:31:52 2017

@author: Niwatori
"""

'''
# Newton's Interpolation method
# Input: Data points, corresponding values and number of points
# Output: Interpolation polynomial
'''
import pylab as pl

def NewtonInterp(x, y, n):
    f = [[0 for i in range(n)] for j in range(n)]
    for i in range(0, n):
        f[i][i] = y[i]

    # Compute divided differences
    for j in range(1, n):
        for i in range(0, n-j):
            f[i][i+j] = (f[i][i+j-1]-f[i+1][i+j]) / (x[i]-x[i+j])
            
    # Compute interpolation polynomial
    def InterpPoly(t):
        ans = f[0][n-1]
        for i in range(n-2, -1, -1):
            ans = ans * (t-x[i]) + f[0][i]
        return ans
    return InterpPoly

# Runge function r(x) and approximate polynomial p(x)
x = range(-5, 6)
r = lambda x:1/(1+x**2)
p = NewtonInterp(x, [r(t) for t in x], 11)

# Plot images of two functions
X = pl.linspace(-5, 5, 100)
pl.plot(X, [r(t) for t in X], label="Runge", linewidth=2)
pl.plot(X, [p(t) for t in X], label="Approximation 1")
pl.legend(loc='upper center')
pl.xlim(-5.5, 5.5)
pl.savefig("fig1.eps")