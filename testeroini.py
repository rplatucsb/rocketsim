# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 23:13:11 2019

@author: Adam
"""
import numpy as np

a = np.array([[1,2],[3,1],[2,2],[3,4]])
b = [a[:,0]==1]
c = [a[:,1]==2]
q = np.logical_and(b,c)
a[q]