#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 22:58:53 2017

@author: Niwatori
"""

''' Spline Interpolation with scipy package '''

import pylab as pl
from scipy.interpolate import interp1d

x = pl.linspace(-5, 5, 11)
r = lambda x:1/(1+x**2)
y = [r(t) for t in x]
     
f = interp1d(x, y, kind = 'slinear')
xnew = pl.linspace(-5, 5, 101)
ynew = f(xnew)
pl.plot(xnew, ynew)