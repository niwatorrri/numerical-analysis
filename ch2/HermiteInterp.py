#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 22:58:53 2017

@author: Niwatori
"""

'''
# Piecewise Cubic Hermite Interpolation
# Input: Data points, corresponding values and derivatives, number of points
# Output: Interpolation function
'''
import pylab as pl

def HermiteInterp(x, y, d, n):
    def InterpPoly(t):
        for i in range(0, n):
            if x[i] <= t <= x[i+1]:
                delta = x[i+1]-x[i]
                ans = y[i] * (1+2*(t-x[i])/delta) * ((t-x[i+1])/delta)**2
                ans += y[i+1] * (1-2*(t-x[i+1])/delta) * ((t-x[i])/delta)**2 
                ans += d[i] * (t-x[i]) * ((t-x[i+1])/delta)**2
                ans += d[i+1] * (t-x[i+1]) * ((t-x[i])/delta)**2
                return ans
    return InterpPoly

x = range(-5, 6)
r = lambda x:1/(1+x**2)
d = lambda x:-2*x/(1+x**2)**2
p = HermiteInterp(x, [r(t) for t in x], [d(t) for t in x], 11)

X = pl.linspace(-5, 5, 199)
pl.plot(X, [r(t) for t in X], label="Runge", linewidth=2)
pl.plot(X, [p(t) for t in X], label="Approximation 4")
pl.legend(loc='upper right')
pl.xlim(-5.5, 5.5)
pl.savefig("fig4.eps")