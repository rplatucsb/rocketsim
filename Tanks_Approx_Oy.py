from math import pi
import math

r1   = 3.25
dens = 0.098
Ysig  = 40000

def tank_para(VCH4,VLOx,P):
    t  = 2*P*r1/(Ysig/1.5)
    r2 = math.sqrt((r1**2)/((VCH4/VLOx)+1))
    L  = VLOx/(pi*(r2**2))
    W  = (pi*L*(((r1+t)**2-r1**2)+((r2+t)**2-r2**2))+(8/3)*pi*((r1+t)**2-r1**2))*dens
    return t,r2,L,W
    
    
