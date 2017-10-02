#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 18:35:26 2017

@author: Niwatori
"""

'''
# Romberg Integration Method
# Input: Target function f, interval [a,b],
#        (Optional) number of subintervals
# Output: Approximate integration of f on [a,b]
'''

def RombergInt(f, a, b, n = 1):
    # Trapezoidal rule
    def T(f, a, b, n):
        ans = 0
        h = (b-a)/n
        for i in range(0, n):
            ans += (f(a+i*h)+f(a+i*h+h))*h/2
        return ans
    
    # Romberg Integration
    M = 15
    R = [[0] * M] * M
    for k in range(1, M):
        R[1][k-1] = T(f, a, b, n*2**(k-1))
        for i in range(1, k):
            R[i+1][k-i-1] = (4**i*R[i][k-i]-R[i][k-i-1])/(4**i-1)
        print('%.16f'%(R[k][0]))
    return R[M-1][0]

# Romberg integration for pi approximation
f = lambda x:4/(1+x*x)
print('%.16f'%(RombergInt(f, 0, 1)))
