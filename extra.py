# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 20:03:20 2024

@author: micha
"""
import numpy as np

def vector_integral(x, y):
    dx = abs(x[0] - x[1])
    return np.sum(y * dx)

def x2(x):
    return x**2

x = np.linspace(0, 10, num=100)
y = x2(x)

print(vector_integral(x, y))
print(np.trapz(y, dx=0.1))
print(0.33 * 10**3)