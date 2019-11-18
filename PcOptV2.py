# -*- coding: utf-8 -*-
"""
Created on Sat May 18 22:27:07 2019

@author: Adam
"""

import numpy as np
import He_mass 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import v5_mars_appropriated as rockSim
import NozzleCreator
import RhaoEq as engineProfile
from math import log10,floor
import Tanks_Approx_W_6_1_19 as propTank

def sigfigs(x,n):
    return round(x, n-int(floor(log10(abs(x)))))

#Get the Specific Impulse interpolations from RPA 
fName = "RPAisp.txt"
f = open(fName,"r")
flines = f.readlines()
flines = flines[8:-1]
rpaPoints = []
for line in flines:
    line = line.replace("\n","")
    elements = line.split(' ')
    elements = list(filter(None,elements))
    data = [float(elements[0]),float(elements[1]),float(elements[10]),float(elements[11])] #MR, Pc, ISP1ATM, ISPVAC
    rpaPoints.append(data)
    
rpaPoints = np.array(rpaPoints)

def interpISP(mr,pc):
    return griddata(rpaPoints[:,0:2],rpaPoints[:,2],(mr,pc),method = 'linear')
    
m = [2.2,
35.0,
3.0,
3.5,
15.0,
1.5,
2.0,
2.0,
5.0,
5.0,
1.0,
4.0,
5.0,
5.0,
0.4,
0.4,
0.5,
0.5,
0.8]

fuelGuesslb = 35
massGuesslb = np.sum(m) - fuelGuesslb
indepMasslb = massGuesslb - ( 1.5 + 2 + 15)

#COM calcs. first num  relative to top of nosecone from spec sheet 6/1/19
comTop,mTop = 2.53,22.7
comTop *= 12 * 2.54
mTop /= 2.2 # in

#define constants
g = 9.81
thrustEst = 3114 #N
lbtoN = 4.44822 #Converts lbf to n
impulse = 9208 #Total Impulse (9208 lb-sec Max)
rhoOx = 1141 #(density of liquid oxygen = 1141 kg/m**3)
rhoMeth = 425.6 #(density of liquid methane = 425.6 kg/m**3
impulseN = impulse * lbtoN

#inputs in psi, dimensionless, lb, in
def apogeeCalc(inputs,independantMass = indepMasslb, rocketDiameter = 6.5,display = False):
    chamberPressure,mixtureRatio = inputs[0],inputs[1]
    
    #pressure losses
    dpInject = .2 * chamberPressure
    dpPipe = .2 * chamberPressure
    tankPressure = chamberPressure + dpInject + dpPipe
    
    parameters,values = [],[]
    
    IspSL = interpISP(mixtureRatio,chamberPressure)
    
    if(display):
        parameters.extend(["ISP at Sea Level (s)"])
        values.extend([IspSL])
    
    #get first order burntime
    tFunc,thrust = NozzleCreator.altThrust(chamberPressure,mixtureRatio)
    apogee,burnTime = rockSim.main(tFunc,rocketDiameter/2,fuelGuesslb,massGuesslb,bTimeCalc = True)
    
    #get engine mass
    abThickness = .909 #cm varies from .906 to .912 from 200 to 500psi at chamber, .923 to .927 at throat, .880 to .885 at end of nozzle
    rhoAb = 1100 # kg/m^3
    xD,yD = engineProfile.design(chamberPressure,mixtureRatio,disp = False,chamberRadiusScaleUpFactor=(1)) #coordinates from RHAO Code
    vIn = np.pi * np.trapz(yD**2,xD) #inner volume
    m = np.array(-1 * np.gradient(yD,xD)) 
    m = np.where(np.isnan(m),0,m)
    xAb,yAb = xD + abThickness/np.sqrt(1+m**2)*m, yD + abThickness/np.sqrt(1+m**2) 
    vAbout = np.pi * np.trapz((yAb)**2,xAb)
    vAb = vAbout - vIn
    vAb *= 1e-6 #cm^3 to m^3
    mAblative = vAb * rhoAb
    eLen = max(xD)-min(xD)
    if(display):
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        plt.plot(xAb,yAb)
        plt.plot(xD,yD)
        ax.set_ylim(-15,15)
        parameters.extend(["Ablative Layer Mass (kg)","Exit Diameter (cm)","Throat Diameter","Engine Length","Chamber Length"])
        values.extend([mAblative,2*max(yD),2*min(yD),eLen,xD[np.argmin(yD)]-xD[0]])
        
    #get propellants masses and volumes
    def propellantData(burnTime,disp = False):
        mDot = thrustEst / (g * IspSL)
        mProps = mDot * burnTime
        mFuel = mProps/(1+mixtureRatio) #fuel mass
        mOxidiser = mProps/(1+1/mixtureRatio) #ox mass
        vOxidiser = mOxidiser/rhoOx #ox volume (density of liquid oxygen = 1141 kg/m**3)
        vFuel = mFuel/rhoMeth #fuel volume (density of liquid methane = 425.6 kg/m**3)
        if(disp):
            parameters.extend(["Mass of Fuel (kg)","Mass of Oxidiser","Mass of Propellants",
                               "Volume of Fuel (l)","Volume of Oxidiser","Volume of Propellants","Mdot (kg/s)"])
            values.extend([mFuel,mOxidiser,mProps,vFuel*1000,vOxidiser*1000,(vFuel+vOxidiser)*1000,mDot])
        return mProps,vOxidiser,vFuel
    
    #get helium tank mass
    def pressureData(disp = False):
        heMass = He_mass.He_mass(vOxidiser,vFuel,tankPressure)
        if(disp):
            parameters.extend(["Helium plus helium tank mass (kg)"])
            values.extend([heMass])
        return heMass
    
    def tankData(vMeth,vLox,disp = False):
        thickness,innerRadius,length,mass = propTank.tank_para(vMeth * 61023.7,vLox * 61023.7,tankPressure)
        if(disp):
            parameters.extend(["Coax tank mass (lb)"])
            values.extend([mass])
        return mass,length
            
    #first run, with guesstimates for masses
    mProps,vOxidiser,vFuel = propellantData(burnTime)
    tankMass,tankLen = tankData(vFuel,vOxidiser)
    heMass = pressureData()
    dryMass = 2.2 * (mAblative + heMass) + independantMass + tankMass
    apogee2,burnTime2 = rockSim.main(tFunc,rocketDiameter/2,mProps * 2.2,dryMass,bTimeCalc = True) 

    #final run, using true-er masses
    mProps,vOxidiser,vFuel = propellantData(burnTime2,disp=display)
    tankMass,tankLen = tankData(vFuel,vOxidiser,disp=display)
    heMass = pressureData(disp = display)
    dryMass = 2.2 * (mAblative + heMass) + independantMass + tankMass
    apogee3,burnTime3 = rockSim.main(tFunc,rocketDiameter/2,mProps * 2.2,dryMass) 
    
    #get center of mass
    #mAblative = 1
    comBot,mBot = 8 * 12 * 2.54 + tankLen + (8.2 *(12*2.54) / 2.2 +(mAblative)*(eLen/2+(2/3*12*2.54)))/(12.3/2.2+(mAblative)),12.3/2.2+(mAblative)
    comNet = (comTop*mTop+comBot*mBot+tankMass*(8*12*2.54+tankLen/2))/(tankMass+mTop+mBot)
    lNet = (3.21 + 2.01 + 2.7)*12*2.54 + tankLen
    
    if(display):
        parameters.extend(["Burn Time (s)","Apogee (ft)","Dry Mass (lb)","Distance Top to Center of mass (cm)","Rocket Length"])
        values.extend([burnTime3,apogee3,dryMass])
    if(display):
        print("       Parameter                   |     Value     ")
        for parameter, value in list(zip(parameters, values)):
            print(parameter + " "*(35-len(parameter))+"| " + str(sigfigs(float(value),3)))
    return apogee3
    
apogeeCalc([350,2.8],display=True)