# -*- coding: utf-8 -*-
"""
Created on Sat May 18 12:49:05 2019

@author: Adam
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

fName = "RPAisp.txt"
f = open(fName,"r")
flines = f.readlines()
flines = flines[8:-1]

rpaPoints = []
mrstep = .05
pcstep = 10
for line in flines:
    line = line.replace("\n","")
    elements = line.split(' ')
    elements = list(filter(None,elements))
    data = [float(elements[0]),float(elements[1]),float(elements[10]),float(elements[11])] #MR, Pc, ISP1ATM, ISPVAC
    rpaPoints.append(data)
    
rpaPoints = np.array(rpaPoints)
MR = np.arange(min(rpaPoints[:,0]),max(rpaPoints[:,0])+mrstep,mrstep)
PC = np.arange(min(rpaPoints[:,1]),max(rpaPoints[:,1])+pcstep,pcstep)
x,y = np.meshgrid(MR,PC)

def interp(mr,pc):
    return griddata(rpaPoints[:,0:2],rpaPoints[:,2],(mr,pc),method = 'linear')


isp = np.zeros((np.size(MR),np.size(PC)))
q = rpaPoints[:,0]
t = rpaPoints[:,1]
for pci in range(np.size(PC)):
    for mri in range(np.size(MR)):
        mrind = np.where(q == np.round(MR[mri],decimals=2))
        pcind = np.where(t == np.round(PC[pci]))
        ind = np.intersect1d(mrind, pcind)
        isp[mri][pci] = rpaPoints[int(ind),2]

MR = np.random.rand(20) + 2.4
PC = np.random.rand(20)*300 + 200
x2,y2 = np.meshgrid(MR,PC)
isp2 = np.zeros((np.size(MR),np.size(PC)))
for pci in range(np.size(PC)):
    for mri in range(np.size(MR)):
        isp2[mri][pci]  = interp(MR[mri],PC[pci]) 

isp = isp.T
isp2 = isp2.T
plt.figure(2)
ax = plt.axes(projection='3d')
ax.set_xlabel("MR)")
ax.set_ylabel("PC")
ax.set_zlabel("ISP")
ax.plot_surface(x,y,isp, cmap='viridis', edgecolor='none')
ax.plot_surface(x2,y2,isp2, cmap='magma', edgecolor='none')

