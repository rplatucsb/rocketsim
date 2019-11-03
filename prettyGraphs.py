# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 20:23:01 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import PcOptV2 as alt
from mpl_toolkits.mplot3d import Axes3D

def f(inputs):
    return -1 * alt.apogeeCalc(inputs)
#print(opt.minimize(f,[250,3],method='trust-constr',bounds = ((150,800),(2.3,3.3)),options={'disp':True,'maxiter' : 5}))
#print(opt.fmin_tnc(f,[200,3],bounds = ((150,800),(2.3,3.3)),approx_grad=True))

count = 10
chamberPressure = np.linspace(150,800,count)
mixtureRatio = np.linspace(2.3,3.6,count)
altitude = np.zeros((count,count,3))

for mri,mr in enumerate(mixtureRatio):
    maxmass = True
    for pci,pc in enumerate(chamberPressure):
        inputs = [pc,mr]
        altitude[mri][pci][0] = pc
        altitude[mri][pci][1] = mr
        altitude[mri][pci][2] = alt.apogeeCalc(inputs,display=False)    
        print(pci+(mri*count))

x,y = np.meshgrid(chamberPressure,mixtureRatio)

plt.figure(2)
ax = plt.axes(projection='3d')
ax.set_xlabel("Chamber pressure (psi)")
ax.set_ylabel("Mixture Ratio ")
ax.set_zlabel("Apogee (ft))")
ax.plot_surface(x,y, altitude[:,:,2], cmap='viridis', edgecolor='none')

t = altitude[np.where(altitude[:,:,2] == np.max(altitude[:,:,2]))]
print(t)



