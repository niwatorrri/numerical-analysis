#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 16:02:48 2017

@author: Niwatori
"""

'''
# Composite Numerical Integration
# Input: Target function f, interval [a,b], number of subintervals n,
#        and the rule for integration (Simpson's rule by default)
# Output: Approximate integration of f on [a,b]
'''
   
def Integrate(f, a, b, n, rule = 'Simpson'):
    ans = 0
    h = (b-a)/n
    for i in range(0, n):
        if rule == 'Midpoint':
            ans += f(a+i*h+h/2)*h
        elif rule == 'Trapezoidal':
            ans += (f(a+i*h)+f(a+i*h+h))*h/2
        elif rule == 'Simpson':
            ans += (f(a+i*h)+4*f(a+i*h+h/2)+f(a+i*h+h))*h/6
        else:
            raise RuntimeError('Invalid rule')
    return ans

# Numerical integration for pi approximation
f = lambda x:4/(1+x*x)
for k in range(1, 7):
    print('$10^{-%d}$&%.16f&%.16f&%.16f\\\\'%(k,
        Integrate(f, 0, 1, 10**(k), rule = 'Midpoint'),
        Integrate(f, 0, 1, 10**(k), rule = 'Trapezoidal'),
        Integrate(f, 0, 1, 10**(k), rule = 'Simpson')))
