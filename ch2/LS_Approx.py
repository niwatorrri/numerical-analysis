#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 19:26:37 2017

@author: Niwatori
"""

'''
# Least Squares Approximation
# Input: Target function, interval and degree of the polynomial
# Output: Approximation polynomial
'''

import pylab as pl
from scipy import zeros, integrate, linalg
from math import sin, pi

def LS_Approx(f, inf, sup, n):
    A = zeros((n+1, n+1))
    base = [(lambda arg:(lambda x:x**arg))(k) for k in range(0, n+1)]
            
    # Compute the equation matrix
    for i in range(0, n+1):
        for j in range(0, n+1):
            A[i][j] = integrate.quad(lambda x:base[i](x)*base[j](x), inf, sup)[0]
    b = [integrate.quad(lambda x:f(x)*base[i](x), inf, sup)[0] for i in range(0, n+1)]
    c = linalg.solve(A, b)
    
    # Compute approximation polynomial
    def ApproxPoly(t):
        ans = 0
        for i in range(0, n+1):
            ans += c[i]*base[i](t)
        return ans
    return ApproxPoly

f = lambda x:sin(pi*x)
p = LS_Approx(f, 0, 1, 2)

X = pl.linspace(0, 1, 99)
pl.plot(X, [f(t) for t in X], label="f(x)=sinπx", linewidth = 2)
pl.plot(X, [p(t) for t in X], label="P*(x)")
pl.legend(loc='upper right')
pl.xlim(-0.1, 1.1)
pl.ylim(-0.1, 1.1)
pl.savefig("fig6.eps")