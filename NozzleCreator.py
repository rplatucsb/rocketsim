# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 18:57:48 2018

@author: Adam

EXIT AREA OFF FROM BOOK BY 4%
"""
import numpy as np
from scipy.interpolate import griddata
#Get the k interpolations from RPA 
fName = "RPAisp.txt"
f = open(fName,"r")
flines = f.readlines()
flines = flines[8:-1]
rpaPoints = []
for line in flines:
    line = line.replace("\n","")
    elements = line.split(' ')
    elements = list(filter(None,elements))
    data = [float(elements[0]),float(elements[1]),float(elements[7])] #MR, Pc, ISP1ATM, ISPVAC
    rpaPoints.append(data)
    
rpaPoints = np.array(rpaPoints)

def kInterp(mr,pc):
    return griddata(rpaPoints[:,0:2],rpaPoints[:,2],(mr,pc),method = 'linear')

#k = 1.1655#dimensionless cnst ratio of specific heats of gas mixture
pOut = .101325  * 10**6 #0.101325 * 10**6 #pressure at exit of nozzle in mpa
#pChamber = 300 * 6894.76 #np.array([100,200,300,400,500]) * (0.00689476)  * 10**6  #first number is in psi
desiredThrust = 3114 #Newtons
thrustCorrectionFactor = .95
g=9.8
ISP = 260 #seconds
pRange = np.linspace(.101325  * 10**6,(.101325/2)*10**6,20)
def dims(pChamber,mRatio,disp):
    k = kInterp(mRatio,pChamber)
    pChamber *= 6894.76
    tCoeff = np.sqrt((2*(k**2))/(k-1)*((2/(k+1))**((k+1)/(k-1)))*(1-(pOut/pChamber)**((k-1)/k))) 
    throatArea = desiredThrust/(pChamber*thrustCorrectionFactor*tCoeff) #m^2 is unit here
    exitAreaRatio = 1/((((k+1)/2)**(1/(k-1)))*((pOut/pChamber)**(1/k))*np.sqrt(((k+1)/(k-1))*(1-((pOut/pChamber)**((k-1)/k)))))
    minChamberArea = 4 * throatArea #book says 3 for appreciable pressure drop, 4 for appreciable chamber velocity
#    mDot = desiredThrust/(ISP*g)
    if(disp):
        print("Thrust coefficient is " , tCoeff)
        print("expansion ratio is ", exitAreaRatio)
        print("throat area is ", throatArea*10000, "cm^2 at tcoeff = ", tCoeff)
        print("exit area is", (throatArea*exitAreaRatio)*10000,"cm^2 and chamber area should be at least ",minChamberArea*10000,"cm^2")
        print(2*np.sqrt(throatArea/np.pi)*100, "cm is the throat diameter", 2*np.sqrt(minChamberArea/np.pi)*100, "cm is the min chamber diameter", 2*np.sqrt(throatArea*exitAreaRatio/np.pi)*100,"cm is exit diameter")
    return exitAreaRatio,2*np.sqrt(throatArea/np.pi)*100,2*np.sqrt(throatArea*exitAreaRatio/np.pi)*100

#give a thrust of altitude curve, acurate within -2 of RPA
def altThrust(pChamber,mRatio):
    k = kInterp(mRatio,pChamber)
    exitArea,throatDiameter,exitDameter = dims(pChamber,mRatio,False)
    pChamber *= 6894.76
    tCoeff = np.sqrt((2*(k**2))/(k-1)*((2/(k+1))**((k+1)/(k-1)))*(1-(pOut/pChamber)**((k-1)/k))) + (pOut-pRange)/pChamber * ((exitDameter/2)**2*np.pi)/((throatDiameter/2)**2*np.pi)
    thrust = tCoeff * ((throatDiameter/200)**2*np.pi) * pChamber * thrustCorrectionFactor
    tTransform = np.polyfit(pRange,thrust,1)
    tFunc = np.poly1d(tTransform)
    return tFunc,thrust

#import matplotlib.pyplot as plt
#tFunc,thrust = altThrust(300,2.8)
#plt.plot(pRange,thrust,'o')
#plt.plot(pRange,tFunc(pRange))
#print(pRange)
#print(thrust)
#dims(300,2.8,True)
#print("AreaRatio", exitAreaRatio)
#print(np.array([600,700,800,900,1000])/14.6959)


