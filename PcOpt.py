 # -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 17:16:22 2019

@author: Adam P, Andrew C
"""
import numpy as np
import v4_mars as rockSim
import RhaoEq as nProf
import scipy.optimize as opt
import tank_solution
import He_mass 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
#inputs in psi,lb,in
def pcOpt(cPress,mass,diameter,mixtureRatio = 2.8, moreData = False):
    nPropMass = mass * 0.453592 #mass of elements of rocket not calculated here
    d = diameter * .0254 #rocket diameter in m
    It = 9208 #Total Impulse (9208 lb-sec Max)
    It = 4.44822*It #Total Impulse (Ns)
    
    parameters,values = [],[]
    #pcD = np.array([200,400,600,800,1000]) data for cnst o/f of 2.8
    #IspD = np.array([2479.3,2736.5,2863.9,2945.6,3004.6])
    #f = np.polyfit(pcD,IspD,4)
    #q = np.poly1d(f) #polynomial least squares fit
    
#    IspD2 = np.array([2169.8,2472.9,2624.9,2721.5,2790.9,2844.2,2887.1,2922.8,2953.2,2979.6]) # ISP CEA data
#    pcD2 = np.array([100,200,300,400,500,600,700,800,900,1000])
#    
#    f2 = np.polyfit(pcD2,IspD2,8)
#    q2 = np.poly1d(f2) #polynomial least squares fit for ISP data
#    
#    
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
    def interpISPV(mr,pc):
        return griddata(rpaPoints[:,0:2],rpaPoints[:,3],(mr,pc),method = 'linear')

    IspSL = interpISP(cPress,mixtureRatio)
    IspVac = interpISPV(cPress,mixtureRatio)
    

    pd_inj = 0.2*pc #pressure drop across injector 
    pd_cool = 0#.2*pc #pressure drop across cool ing channels
    pd_plumb = 0.2*pc #pressure drop across plumbing
    pt = pc+pd_inj+pd_cool+pd_plumb #pressure in fuel and ox tanks
    pt_si = pt*6894.757 #chamber pressure (Pa)
    pr = pt_si * 1.5 #rated pressure (Pa) --> safety factor of 1.5

    def propellantsMassVols(disp = False):
        mp = (It)/Isp #propellant mass
        mf = mp/(1+mixtureRatio) #fuel mass
        mo = mp/(1+1/mixtureRatio) #ox mass
        Vo = mo/1141 #ox volume (density of liquid oxygen = 1141 kg/m**3)
        Vf = mf/425.6 #fuel volume (density of liquid methane = 425.6 kg/m**3)
        mp = mo+mf #propellant mass
        if(disp):
            parameters.extend(["Mass of Fuel (kg)","Mass of Oxidiser","Mass of Propellants",
                               "Volume of Fuel (l)","Volume of Oxidiser","Volume of Propellants","Mdot (kg/s)"])
            values.extend([mf,mo,mp,Vf*1000,Vo*1000,(Vf+Vo)*1000,mp/(9208/700)])
        return mf,mo,Vo,Vf
    
#    def propellantTanksMass(oxVol, methVol, disp = False):
#        th_t = (pr*(d/2))/(s_al) #thickness of prop/ox tanks required
#        
#        A_t = np.pi*((d/2)-th_t)**2 #inner cross sectional area for fuel and ox
#        l_ot = oxVol/A_t #length of ox tank
#        l_ft = methVol/A_t #length of fuel tank
#        
#        V_ot = 2*(np.pi*(d/2)**2)*th_t + (np.pi*((d/2)**2-(d/2-th_t)**2))*l_ot #material volume of ox tank
#        V_ft = 2*(np.pi*(d/2)**2)*th_t + (np.pi*((d/2)**2-(d/2-th_t)**2))*l_ft #material volume of fuel tank
#        
#        mot = V_ot*rho_al #ox tank mass
#        mft = V_ft*rho_al #fuel tank mass
#        mt = mot+mft #total tank mass
#        if(disp):
#            parameters.extend(["Oxidiser Tank Mass (kg)","Fuel Tank Mass","Total Tank Mass","Tank Thickness(cm)"])
#            values.extend([mot,mft,mt,th_t*100])
#        return mt
#    
#    def pressurizingSysMass(oxVol,methVol,disp = False):
#        mtHe = []
#        mHe = []
#        vOut = (4/3)*np.pi*(d/2)**3
#        def volumeHe(pressHe):
#            return (4/3)*np.pi*((d/2)*(1-pressHe*1.5/(2*s_al)))**3
#        def hePressure(press,ind):
#            return pr[ind]*((oxVol[ind]+methVol[ind])+volumeHe(press))-press*volumeHe(press)
#        for i in range(0,len(pc)):
#            pHe = opt.ridder(hePressure,10000,100000000,args=i)
#            vHe = volumeHe(pHe)
#            mHe.append((vHe*pHe/(rHe*300))*4)
#            mtHe.append((vOut-vHe) * rho_al)
#        if(disp):
#            parameters.extend(["He Volume(l)","He Pressure (MPA)","He Tank Mass(kg)","He Mass"])
#            values.extend([vHe,pHe,mtHe,mHe])
#        return mtHe + mHe
    
    rhoAb = 2000 #density(kg/m^3) estimate for carbon layer 
    abThickness = 0.8 #thickness (cm) estimate, 
    #above data from https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-512-rocket-propulsion-fall-2005/lecture-notes/lecture_10.pdf
    def engineMass(disp = False):
        mEng = []
        tEng = []
        for i in pc:
            xD,yD = nProf.design(i,disp = disp) #coordinates from RHAO Code
            innerSurface = 2 * np.pi * np.trapz(yD,xD)
            xD += abThickness
            yD += abThickness
            vIn = np.pi * np.trapz(yD**2,xD) #inner volume
            engineThickness = (i*6894.757 * 1.5 * max(yD) * 2 )/s_cu
            m = np.array(-1 * np.gradient(yD,xD)) 
            m = np.where(np.isnan(m),0,m)
            xOut,yOut = xD + engineThickness/np.sqrt(1+m**2)*m, yD+engineThickness/np.sqrt(1+m**2)# outside wall coords
            xAb,yAb = xD - abThickness/np.sqrt(1+m**2)*m, yD-abThickness/np.sqrt(1+m**2) 
            vOut = np.pi * np.trapz((yOut)**2,xOut) #outside volume
            vAbi = np.pi * np.trapz((yAb)**2,xAb)
            vAb = vIn - vAbi
            vAb *= 1e-6
            mAb = vAb * rhoAb
            vEng = vOut - vIn 
            vEng *= 1e-6
            mEng.append(vEng * rho_cu)
            tEng.append(engineThickness)
        if(disp):
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            plt.plot(xAb,yAb)
            plt.plot(xOut,yOut)
            ax.set_ylim(-15,15)
            parameters.extend(["Engine Mass (kg)","Engine Thickness (cm)","Potential Ablative Layer Mass (kg)","Exit Diameter (cm)","Throat Diameter","Engine Length","Chamber Length"])
            values.extend([mEng,tEng,mAb,2*max(yD),2*min(yD),max(xD)-min(xD),xD[np.argmin(yD)]-xD[0]])
        return mEng + mAb
    
    apogee = []
    methM,oxM,oxVol,methVol = propellantsMassVols(disp=moreData)
    mHe = He_mass.He_mass(oxVol,methVol,pr)
    mEng = [0] #engineMass(disp=moreData)
    coaxTank = tank_solution.TankSolution(pr[0]/6894.757,methM[0],oxM[0],disp=moreData)
    #print(propellantTanksMass(oxVol,methVol))
    mt = 0 #(coaxTank.mass_tank_lox + coaxTank.mass_tank_meth)
    if(moreData):
        parameters.extend(["Lox Tank Mass (kg)","Meth Tank Mass"])
        values.extend([coaxTank.mass_tank_lox,coaxTank.mass_tank_meth])
    for i in range(0,len(pc)):
        apogee.append(rockSim.main(700,diameter/2,(methM[0] + oxM[0])*2.2,(nPropMass+mt+mHe)*2.2))
#    print((methM[0] + oxM[0])*2.2,(nPropMass+mt+mEng[0])*2.2)
    dryMass =  nPropMass+mt+mEng[0] + mHe
    wetMass = dryMass + methM[0] + oxM[0]
    if(not moreData):
        return apogee[0]
    else:
        parameters.extend(["Apogee (ft)","Dry Mass(lb)","Wet Mass","Isp (s)"])
        values.extend([apogee,dryMass*2.2,wetMass*2.2,Isp/9.8])
        print("       Parameter                   |     Value     ")
        for parameter, value in list(zip(parameters, values)):
            print(parameter + " "*(35-len(parameter))+"| " + str(value))
        return apogee,(methM+oxM)*2.2 #rockSim.main(700,diameter/2,(methM[0] + oxM[0])*2.2,(nPropMass+mt+mEng[0])*2.2)    
   #,dryMass * 2.2,wetMass*2.2,Isp/9.8

#pc = np.arange(200,400,10.)
#alt = []
#for i in pc:
#    alt.append(pcOpt(i,40,6.5,moreData=False))
#plt.plot(pc,alt)

print(pcOpt(300,5+2.2+3.5+2+2+2+12+2+1.5+2+4+6+3+2,6.5,moreData=True))
#print(v1.main(3114,3.25,16.3,25))   
#print("Max apogee ", max(apogee), " ft at Pc ",  pc[np.argmax(apogee)], " psi dry mass ", dryMass , " kg wet mass ", wetMass)
#plt.plot(pc,apogee)
#plt.ylabel("Apogee (ft)")
#plt.xlabel("Pc (psi)")

