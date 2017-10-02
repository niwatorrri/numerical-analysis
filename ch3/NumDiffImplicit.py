#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 00:46:43 2017

@author: Niwatori
"""

'''
# Numerical Differentiation (Implicit form)
# Input: Target function f, desired x's, number of x's and stride h
# Output: Approximate derivatives of f at each x
'''

from math import exp

def DiffExplicit(f, x, h):
    return (f(x+h)-f(x-h))/(2*h)
    
def DiffImplicit(f, x, n, h):
    a = [1]*(n-1)+[0]
    b = [1]+[4]*(n-1)+[1]
    c = [0]+[1]*(n-1)
    d = [DiffExplicit(f, x[0], h)] \
         + [3*(f(x[i+1])-f(x[i-1]))/h for i in range(1, n)] \
         + [DiffExplicit(f, x[n], h)]
    return SolveTridiag(a, b, c, d, n+1)

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

# Comparison of explicit form and implicit form    
f = lambda x:exp(-x/4)
h = 0.1
n = 10
x = [i*h for i in range(0, n+1)]
ex = [DiffExplicit(f, t, h) for t in x]
im = DiffImplicit(f, x, n, h)
ac = [-f(t)/4 for t in x]
for i in range(1, n):
    print('%.1f&%.8f&%.8f&%.8f&'%(i*h, ex[i], im[i], ac[i]))
    print('%.3e&%.3e\\\\'%(abs(ex[i]-ac[i]), abs(im[i]-ac[i])))
