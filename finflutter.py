import numpy as np
import matplotlib.pyplot as plt
import math

G=1.45e6; #Effective Shear Modulus
Altitude=np.linspace(1,45000,45000)
Length1=24
Length2=0
Thickness=.3
Base=8.4
S=(Length1+Length2)/2*Base
AR=Base**2/S
Lam=Length2/Length1
T = []
P = []
#altitude below 36152
for i in range(0,36152,1):
    Ti=59-.00356*Altitude[i]
    T.append(Ti)
    Pi=(2116/144)*((Ti+459.7)/518.6)**5.256
    P.append(Pi)

for i in range(36152,45000,1):
    Ti=-70
    T.append(Ti)
    Pi=473.1/144*math.exp(1.73-.000048*Altitude[i])
    P.append(Pi)

a = []
V = []
M = []

for i in range(0,len(T),1):
    A=math.sqrt(1.4*1716.59*(T[i]+460))
    v=A*math.sqrt((G)/((1.337*AR**3*P[i]*(Lam+1))/(2*(AR+2)*(Thickness/Length1)**3)))
    m = v/A
    a.append(A)
    V.append(v)
    M.append(m)

# plt.plot(Altitude,M)
# plt.xlabel('Altitude')
# plt.ylabel('Mach Flutter')
# plt.show()
