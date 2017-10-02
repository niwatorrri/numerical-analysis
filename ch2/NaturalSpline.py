#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 22:58:53 2017

@author: Niwatori
"""

'''
# Natural Cubic Spline Interpolation
# Input: Data points, corresponding values and number of points
# Output: Interpolation function
'''
import pylab as pl

def SolveTridiag(a, b, c, d, n):
    '''
    [ b c       ] [ x ]   [ d ]
    [ a b c     ] [ x ]   [ d ]
    [   . . .   ] [ x ] = [ d ]
    [     a b c ] [ x ]   [ d ]
    [       a b ] [ x ]   [ d ]
    '''
    a = [0] + a
    for i in range(1, n):
        b[i] -= a[i]*c[i-1]/b[i-1]
        d[i] -= a[i]*d[i-1]/b[i-1]

    x = [0] * n    
    x[n-1] = d[n-1]/b[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (d[i]-c[i]*x[i+1])/b[i]
    return x

def NaturalSpline(x, y, n):
    n -= 1
    delta = [1/(x[1]-x[0])]
    lam = [1]
    mu = [3*delta[0]*(y[1]-y[0])]
    
    # Compute the equation matrix
    for i in range(1, n):
        delta.append(1/(x[i+1]-x[i]))
        lam.append(delta[i]/(delta[i-1]+delta[i]))
        mu.append(3*(1-lam[i])*delta[i-1]*(y[i]-y[i-1])+3*lam[i]*delta[i]*(y[i+1]-y[i]))
    mu.append(3*delta[n-1]*(y[n]-y[n-1]))
    _lam = [1-t for t in lam[1:n]] + [1]
    
    # Solve derivatives at each point
    m = SolveTridiag(_lam, [2] * (n+1), lam, mu, n+1)
    
    # Compute the natural spline function
    def InterpPoly(t):
        for i in range(0, n):
            if x[i] <= t <= x[i+1]:
                d = x[i+1]-x[i]
                ans = y[i] * (1+2*(t-x[i])/d) * ((t-x[i+1])/d)**2
                ans += y[i+1] * (1-2*(t-x[i+1])/d) * ((t-x[i])/d)**2 
                ans += m[i] * (t-x[i]) * ((t-x[i+1])/d)**2
                ans += m[i+1] * (t-x[i+1]) * ((t-x[i])/d)**2
                return ans
    return InterpPoly

x = range(-5, 6)
r = lambda x:1/(1+x**2)
p = NaturalSpline(x, [r(t) for t in x], 11)

X = pl.linspace(-5, 5, 199)
pl.plot(X, [r(t) for t in X], label="Runge function", linewidth = 2)
pl.plot(X, [p(t) for t in X], label="Natural spline")
pl.legend(loc='upper right')
pl.xlim(-5.5, 5.5)
pl.savefig("fig5.eps")