#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:51:56 2017

@author: Niwatori
"""

'''
# Adaptive Quadrature Method
# Input: Target function f, interval [a,b], and the tolerance eps
# Output: Approximate integration of f on [a,b]
'''

def I(f, a, b):
    return (f(a)+4*f((a+b)/2)+f(b))*(b-a)/6

def Integrate(f, a, b, eps, lev = 1):
    tmp = I(f, a, (a+b)/2) + I(f, (a+b)/2, b)
    if (abs(I(f, a, b)-tmp) < 10*eps) or (lev > 15):
        return tmp
    else:
        return Integrate(f, a, (a+b)/2, eps/2, lev+1)\
                + Integrate(f, (a+b)/2, b, eps/2, lev+1)

# Adaptive quadrature for pi approximation
f = lambda x:4/(1+x*x)
for k in range(1, 24, 3):
    print('$10^{-%d}$&%.20f\\\\'%(k, Integrate(f, 0, 1, 10**(-k))))
