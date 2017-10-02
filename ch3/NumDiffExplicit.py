#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 19:26:37 2017

@author: Niwatori
"""

'''
# Numerical Differentiation (Explicit form)
# f1, f2 and f3 are all feasible approximations
# Input: Target function f, desired x and stride h
# Output: Approximate derivative of f at x
'''

from math import log

def f1(f, x, h):
    return (f(x+h)-f(x))/h
def f2(f, x, h):
    return (f(x+h)-f(x-h))/(2*h)
def f3(f, x, h):
    return (f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/(12*h)

f = lambda x:log(x)
x = 0.7
h = 1
ans = 10/7
for k in range(1, 11):
    h *= 0.1
    print('$10^{-%d}$&%.12f&%.12f&%.12f&'%
          (k, f1(f,x,h), f2(f,x,h), f3(f,x,h)))
    print('%.4e&%.4e&%.4e\\\\'%
          (abs(ans-f1(f,x,h)), abs(ans-f2(f,x,h)), abs(ans-f3(f,x,h))))
