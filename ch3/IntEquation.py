#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 14:35:48 2017

@author: Niwatori
"""

'''
# Solve an integral equation using composite Simpson's rule
# \int_0^1 (s^2+t^2)^{1/2}u(t)dt = ((s^2+1)^{3/2}-s^3)/3
# Output: Solution of order n and condition number of matrix
'''

from numpy import zeros, linalg

def Solve(n):
    A = zeros((2*n+1, 2*n+1))
    b = zeros((2*n+1, 1))
    r = lambda i,j:(i*i+j*j)**0.5/(2*n)
    for i in range(0, 2*n+1):
        for j in range(0, 2*n+1):
            A[i][j] = r(i, j)
            if 1 <= j <= 2*n-1:
                if j % 2 == 1: A[i][j] *= 4
                else: A[i][j] *= 2
        b[i] = ((1+(i/(2*n))**2)**1.5-(i/(2*n))**3)*(2*n)
    return [linalg.solve(A, b), linalg.cond(A)]

for k in range(7, 10):
    sol, cond = Solve(k)
    for i in range(0, 2*k+1):
        if i == k: print(k)
        print('&%.8f&%.8f&'%(sol[i], i/(2*k)),  end='')
        if i == k: print('%.4e'%(cond))
        print('\\\\')
    print('\\hline')