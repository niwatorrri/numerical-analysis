#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 17:12:46 2017

@author: Niwatori
"""

'''
# Linear spline interpolation
# Input: Data points, corresponding values and number of points
# Output: Interpolation function
'''
import pylab as pl

def LinearSpline(x, y, n):
    def InterpPoly(t):
        for i in range(0, n):
            if x[i] <= t <= x[i+1]:
                return (y[i+1]*(t-x[i]) - y[i]*(t-x[i+1])) / (x[i+1]-x[i])
    return InterpPoly

x = range(-5, 6)
r = lambda x:1/(1+x**2)
p = LinearSpline(x, [r(t) for t in x], 11)

X = pl.linspace(-5, 5, 199)
pl.plot(X, [r(t) for t in X], label="Runge", linewidth=2)
pl.plot(X, [p(t) for t in X], label="Approximation 3")
pl.legend(loc='upper right')
pl.xlim(-5.5, 5.5)
pl.savefig("fig3.eps")