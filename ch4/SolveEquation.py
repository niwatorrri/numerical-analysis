#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:46:08 2017

@author: Niwatori
"""

'''
# Solve Nonlinear Equations
# Bisection method & Fixed-point iteration (Newton's method included)
# Output: Approximate root
'''

import pylab as pl
from math import sin, cos, acos, exp, log, copysign
sgn = lambda x: copysign(1, x)

''' Bisection Method '''
def Bisect(f, a, b, eps = 1e-8):
    fa = f(a)
    fb = f(b)
    if sgn(fa) * sgn(fb) > 0:
        return 'f(a)f(b)<0 not satisfied.'

    while (b - a) / 2 > eps:
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) < 1e-8:
            break
        if sgn(fa) * sgn(fc) < 0:
            b = c; fb = fc
        else:
            a = c; fa = fc
    return (a + b) / 2

''' Fixed-Point Iteration (or Newton's method) '''
def FP_Iteration(f, x, eps = 1e-8, steps = 300):
    x1 = f(x)
    while abs(x1-x) > eps and steps:
        # print('%d&%.12f&%.4e\\\\'%(301-steps, x1, abs(x1-3.076421163792758)))
        x = x1; x1 = f(x); steps -= 1
    if steps == 0:
        return 'Iteration diverges.'
    return x1

# Solve zeros of f(x)=sin(10x)-x
f = lambda x:sin(10*x)-x
print(Bisect(f, 0, 0.5))
print(Bisect(f, 0.5, 0.8))
print(Bisect(f, 0.8, 1))

g = lambda x:x-(sin(10*x)-x)/(10*cos(10*x)-1)
print(FP_Iteration(g, 0.3))
print(FP_Iteration(g, 0.7))
print(FP_Iteration(g, 0.8))

# Solve root of cos(x)+1/(1+exp(-2x))=0
p1 = lambda x:acos(-1/(1+exp(-2*x)))
p2 = lambda x:0.5*log(-1/(1+1/cos(x)))
p3 = lambda x:x-(cos(x)+1/(1+exp(-2*x)))/(-sin(x)+2*exp(-2*x)/(1+exp(-2*x))**2)
print(FP_Iteration(p1, 3, 1e-12))
# print(FP_Iteration(p2, 3, 1e-12)) # Diverges
print(FP_Iteration(p3, 3, 1e-12))

# Plot images of the functions
X = pl.linspace(-2, 2, 99)
pl.plot(X, [sin(10*t) for t in X], label="y=sin(10x)")
pl.plot(X, [t for t in X], label="y=x")
pl.legend(loc='lower right')
pl.show()
# pl.savefig("fig1.eps")