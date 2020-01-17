# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:33:31 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
import NozzleCreator as nozz
import csv

intocm = 0.393701
#####Data and Rao process from http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
#inputs
#    bevelScale  #length fraction of chamber which will become bezeir curve
#    minT  # angle at which the first arc starts
#    res #resolution
def design(pC,MR,disp = False,res = .001,chamberRadiusScaleUpFactor = 1,bevelScale = .4,minT = -135):
    #Data 
    eRat = np.array([3.5,4,5,6,7,8,9,10,20])
    tNData = np.array([19.8,21.6,23,24,24.8,25.3,26,27.2,28.9])
    tEData = np.array([14.8,14,13,12,11.4,11,10.8,10.2,9])
    tNC = np.polyfit(eRat,tNData,2)
    tNF = np.poly1d(tNC)
    tEC = np.polyfit(eRat,tEData,2)
    tEF = np.poly1d(tEC)
    lStar = 100 #cm
    #nozzle and throat characteristics from my code
    epsilon,dT,dE = nozz.dims(pC,MR,disp)
    rE,rT = dE/2,dT/2
    lN = 0.8*(rT*(np.sqrt(epsilon)-1))/np.tan(np.deg2rad(15))
    tN,tE = tNF(epsilon),tEF(epsilon)
    #all inputs in degrees and cm, outputs in cm
    #define bezeir converging curve
    def conv(t):
        eX,eY = entranceT(minT)
        nX,nY = xC[0] + lCh * (1-bevelScale), rE * chamberRadiusScaleUpFactor
        qY =  rE * chamberRadiusScaleUpFactor
        qX =  eX-(rE * chamberRadiusScaleUpFactor - eY)/(np.tan(np.deg2rad(-minT-90)))
        plt.clf()
        plt.figure(2)
        plt.plot([nX,qX,eX],[nY,qY,eY])
        x = (1-t)**2 * nX + 2*(1-t)*t*qX + t**2*eX
        y = (1-t)**2 * nY + 2*(1-t)*t*qY + t**2*eY  
        return x,y
    
    #Define first arc
    def entranceT(theta): #converging
        theta = np.deg2rad(theta)
        x = 1.5 * rT * np.cos(theta) 
        y = 1.5 * rT * np.sin(theta) + 2.5 * rT
        return x,y
    
    #Define second  arc
    def exitT(theta): #diverging
        theta = np.deg2rad(theta)
        x = .382 * rT * np.cos(theta)
        y = .382 * rT * np.sin(theta) + 1.382 * rT
        return x,y
    
    #Define bezier diverging portion
    def nozzle(t): #nozzle
        nX,nY = exitT(tN-90)
        eX,eY = lN,rE
        m1,m2 = np.tan(np.deg2rad(tN)),np.tan(np.deg2rad(tE))
        c1,c2 = nY - m1*nX, eY - m2*eX
        qX,qY = (c2-c1)/(m1-m2),(m1*c2 - m2*c1)/(m1-m2)
        x = (1-t)**2 * nX + 2*(1-t)*t*qX + t**2*eX
        y = (1-t)**2 * nY + 2*(1-t)*t*qY + t**2*eY  
        return x,y
    
    #find theta for which converging width = nozzle exit width
#    def diff(theta):
#        x1,y1 = entranceT(theta)
#        x2,y2 = nozzle(1)[0],chamberRadiusScaleUpFactor*rE
#        return y1-y2
#       
#    for i in np.arange(-300,-100,.5):
#        if diff(i) < 0 :
#            minT =  i
#            break
    thetaE = np.arange(minT,-90+res,res)
    x,y = entranceT(thetaE)
    
    thetaEx = np.arange(-90, tN - 90+res, res)
    x1,y1 = exitT(thetaEx)
    
    t = np.arange(0,1+res*10,res*10)
    x2,y2 = nozzle(t)
    
    vCon = np.pi * np.trapz(y**2,x)
    vCyl = lStar*(rT**2*np.pi) - vCon
    lCh = vCyl/((chamberRadiusScaleUpFactor*rE)**2*np.pi)
    xC,yC = [x[0]-lCh],[rE * chamberRadiusScaleUpFactor]
    
    xB,yB = conv(t)
    
    repeats = 15
    for i in range(repeats):
        vCon = np.pi * np.trapz(y**2,x)
        vCon += np.pi * np.trapz(yB**2,xB)
        vCyl = lStar*(rT**2*np.pi) - vCon
        lCh = vCyl/((chamberRadiusScaleUpFactor*rE)**2*np.pi)
        xC,yC = np.array([xB[0]-lCh]),np.array([rE * chamberRadiusScaleUpFactor])
        xB,yB = conv(t)

    xNet = np.concatenate((xC,xB,x,x1,x2))
    yNet = np.concatenate((yC,yB,y,y1,y2))
    
    xInterp = np.arange(xC[0],x2[-1],res)
    yInterp = np.interp(xInterp,xNet,yNet)

    def plot():
        plt.figure(1)
        plt.plot(x,y)
        plt.plot(x1,y1)
        plt.plot(x2,y2)
        plt.plot(xB,yB)
        plt.xlabel("X (cm)")
        plt.ylabel("Y (cm)")
        plt.ylim((0,min(yNet) + max(xNet)-min(xNet)))
        plt.figure(2)
        plt.plot(xB,yB)
        plt.plot(x,y)
        plt.figure(3)
        plt.plot(xNet*intocm,yNet*intocm)
        plt.plot(xNet*intocm,-yNet*intocm)
        plt.xlabel("X (in)")
        plt.ylabel("Y (in)")
        plt.ylim(-(max(xNet)-min(xNet))*intocm/2,(max(xNet)-min(xNet))*intocm/2)

    if(disp):
       plot()
       print("Nozzle Length " + str (lN) + " \nTheta N " +  str(tN) + " \nTheta E Chamber " + str(tE) + " \nCylinder length/diameter " + str(lCh/(2 * chamberRadiusScaleUpFactor*rE)) + " \nchamber length " + str(max(x)-min(xC)))
       print("Chamber L/D " + str(-xInterp[0]/(2*rE*chamberRadiusScaleUpFactor)))
       pass
    return xInterp,yInterp


xVals,yVals = design(350,2.8,disp=True,chamberRadiusScaleUpFactor=1.65,bevelScale=.04)
yC = np.where(yVals == min(yVals))
cVol = np.pi * np.trapz(yVals[:int(yC[0])]**2,xVals[:int(yC[0])])
print("Chamber L* : " + str(cVol/(min(yVals)**2*np.pi)))
print("Contraction Ratio : " + str(max(yVals)**2/min(yVals)**2))
#Code for Generating A csv or txt of our nozzle
#xVals *= .01
#xVals -= xVals[0]
#yVals *= .01
#
#file = open(r'nozzfile.txt','w')
#for i in range(len(xVals)):
#    file.write(str(xVals[i]) + " " + str(yVals[i])+' 0 \n')
#file.close()
#with open('Nozzle350psi2.8MR.csv',mode='w',newline='') as nozzFile:
#    writer = csv.writer(nozzFile)
#    writer.writerow(xVals)
#    writer.writerow(yVals)
#
