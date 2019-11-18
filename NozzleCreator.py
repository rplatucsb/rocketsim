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

pOut = .101325  * 10**6 #0.101325 * 10**6 #pressure at exit of nozzle in mpa
desiredThrust = 3114 #Newtons
thrustCorrectionFactor = .9527 # estimated efficiency of our nozzle
psiToMpa = 6894.76 
#Find engine areas
def dims(pChamber,mRatio,disp=False):
    k = kInterp(mRatio,pChamber)
    pChamber *= psiToMpa
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
    pRange = np.linspace(.101325  * 10**6,(.101325/2)*10**6,20)
    k = kInterp(mRatio,pChamber)
    exitArea,throatDiameter,exitDameter = dims(pChamber,mRatio,False)
    pChamber *= psiToMpa
    tCoeff = np.sqrt((2*(k**2))/(k-1)*((2/(k+1))**((k+1)/(k-1)))*(1-(pOut/pChamber)**((k-1)/k))) + (pOut-pRange)/pChamber * ((exitDameter/2)**2*np.pi)/((throatDiameter/2)**2*np.pi)
    thrust = tCoeff * ((throatDiameter/200)**2*np.pi) * pChamber * thrustCorrectionFactor
    tTransform = np.polyfit(pRange,thrust,1)
    tFunc = np.poly1d(tTransform)
    return tFunc,thrust

