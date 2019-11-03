# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:58:44 2019

@author: Adam
"""
import numpy as np
a = 1262 #m/s
D = .07165 
L = .247

k = 0  #long
Rad = 0
Tang = 1
lambd =  0
if(Rad == 1 and Tang == 0):
    lambd = 3.8317
if(Tang == 1 and Rad == 0):
    lambd = 1.8412
if(Tang == 2 and Rad == 0):
    lambd = 3.0541
f2 = a/(2*np.pi) * np.sqrt(lambd**2/((D/2)**2)+(k*np.pi/L)**2)
print(f2)