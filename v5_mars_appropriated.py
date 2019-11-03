"""

UCSB Rocket Propulsion Laboratory
Flight Dynamics Rocket Simulation

Andrew Zakoor
Nolan McCarthy
Adam Poklemba

Last Edit: 05/07/2019 3:06PM ANDREW ZAKOOR

"""

from math import pi,pow,exp,sqrt,sin,cos,tan,atan
from matplotlib.font_manager import FontProperties
from sympy.solvers import solve
from sympy import Symbol
import threading
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import random
import sys

args = sys.argv

def main(thrustofAlt,innerRadius_in,fuel_mass_lb,dry_rocket_mass_lb,bTimeCalc = False):

    def rotation(wind_speed,theta):

        # first run through, initial data

        i6_windspeeds = [0, 75, 100, 150, 200, 225, 300, 375, 400, 450, 525, 600] #crosswinds in in/s (corresponds to 0 --> 30 mph)
        i6_windtorque = [0, 2.6781, 4.9022, 11.378, 20.4606, 26.1173, 47.3316, 62, 75.2607, 85.1719, 109.523, 148.944, 195.733] #torques into wind in lbf*in

        i6_headwinds = [2040, 11500] #in/s, absolute value of 10 degree a.o.a.
        i6_restorque = [688.526, 8745.068] #lbf*in, around CG at 59.39 in from bottom
        i6_windforce = [25.806, 344.834] #CP is found by torque/force

        # second run through, refined, more accurate data

        i6_headtorque_rail_bait = [-693.025, -508.286, -323.683, -184.562, -130.647, -92.529, -35.430, 0, 35.430, 92.529, 130.647, 184.562, 323.683, 508.286, 693.025] #restoring torque from headwinds (lbf*in)
        i6_headtorque_rail =[i*((velocity)/(48.78))**1.4 for i in i6_headtorque_rail_bait] #multiplies torque by ratio
        i6_headangles_rail_bait = [-10, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10] #corresponding angles assumed for rocket coming off rails
        i6_headangles_rail = [i*(pi/180) for i in i6_headangles_rail_bait] #in radians

        i6_windtorque_30_bait = [137.537, 142.872, 147.621, 150.783, 150.288, 150.783, 147.621, 142.872, 137.537] #displacing torque from wind, for speeds of 30mph at varying degrees (lbf*in)
        i6_windtorque_30 = [i*(wind_speed/13.41) for i in i6_windtorque_30_bait] #converts to ratio
        i6_windangles_30_bait = [-10, -6, -3, -1, 0, 1, 3, 6, 10] #angles for rocket at 30mph winds
        i6_windangles_30 = [i*(pi/180) for i in i6_windangles_30_bait] #in radians

        # I_rocket = (1/12)*(fuel_mass_lb+dry_rocket_mass_lb)*(rocket_length)**2 + (dry_rocket_mass+fuel_mass_lb)*(cg_list[0]-rocket_length/2)**2
        I_rocket = 131014.23 #lbs*in^2

        # while omega <= 0:

        for i in range(0,len(i6_headangles_rail),1):
            if theta >= i6_headangles_rail[i] and theta < i6_headangles_rail[i+1]:
                def coefficient(x,y):
                    x_1 = x[0]
                    x_2 = x[1]
                    x_3 = x[2]
                    y_1 = y[0]
                    y_2 = y[1]
                    y_3 = y[2]

                    a = y_1/((x_1-x_2)*(x_1-x_3)) + y_2/((x_2-x_1)*(x_2-x_3)) + y_3/((x_3-x_1)*(x_3-x_2))

                    b = (-y_1*(x_2+x_3)/((x_1-x_2)*(x_1-x_3))
                         -y_2*(x_1+x_3)/((x_2-x_1)*(x_2-x_3))
                         -y_3*(x_1+x_2)/((x_3-x_1)*(x_3-x_2)))

                    c = (y_1*x_2*x_3/((x_1-x_2)*(x_1-x_3))
                        +y_2*x_1*x_3/((x_2-x_1)*(x_2-x_3))
                        +y_3*x_1*x_2/((x_3-x_1)*(x_3-x_2)))

                    return a,b,c

                x = [i6_headangles_rail[i],i6_headangles_rail[i+1],i6_headangles_rail[i+2]]
                y = [i6_headtorque_rail[i],i6_headtorque_rail[i+1],i6_headtorque_rail[i+2]]

                a,b,c = coefficient(x,y)

                head_torque_q = a*(theta**2)+b*(theta)+c
                break

        for j in range(0,len(i6_windangles_30),1):
            if theta >= i6_windangles_30[j] and theta < i6_windangles_30[j+1]:
                m4 = (i6_windtorque_30[j+1]-i6_windtorque_30[j])/(i6_windangles_30[j+1]-i6_windangles_30[j])
                wind_torque = (theta-i6_windangles_30[j]) * m4 + i6_windtorque_30[j]
                break

        current_torque = head_torque_q - wind_torque
        alpha = current_torque/I_rocket

        return theta, omega, alpha

    def cd_altitude_calc():

        slope_a = ( 224 - 293 ) / 11000 #kelvin decrease per meter
        a_slope = 0.0065 #the slope of temperature decrease K/m for h<11000 meters

        if height < 11000: #calculates the temperature at a given height
            altitude_temperature = mojave_temperature - ( height * meters_to_feet * temp_slope )
            altitude_density = mojave_density * (mojave_temperature/(mojave_temperature-a_slope*height))**(1+(avg_g*M)/(-a_slope*gas_R))
            altitude_pressure = mojave_pressure * (mojave_temperature/(mojave_temperature-a_slope*height))**(avg_g*M/(-a_slope*gas_R))
        else: #above 11k meters the temperature is relatively constant
            altitude_temperature = mojave_temperature - ( 11000 * meters_to_feet * temp_slope )
            altitude_density = mojave_density * math.exp((-avg_g*M*height)/(gas_R*mojave_temperature))
            altitude_pressure = mojave_pressure * math.exp((-avg_g*M*height)/(gas_R*mojave_temperature))

        altitude_mach = sqrt( gamma * gas_constant * altitude_temperature ) #mach from temperature, in m/s

        cd_v3_list = np.array([0.55432, 0.53437, 0.51071, 0.48687, 0.49139, 0.48243, 0.47812, 0.47577, 0.46904, 0.49098, 0.51463, 0.68300, 0.67543, 0.62601, 0.52901, 0.46041, 0.37132, 0.35000])
        drag_v3_list = np.array([0, 8.77484, 33.31318, 72.91437, 127.70369, 197.82749, 283.87060, 387.60258, 511.45387, 579.75506, 691.77405, 1084.06154, 1267.29417, 1426.57258, 1969.45133, 2425.65958, 3757.03066])
        drag_v4_list = np.array([0, 8.52192, 32.28034, 70.80437, 123.82969, 191.56966, 274.50500, 374.45160, 493.56528, 671.10050, 833.52275, 1338.57703, 1416.65956, 1556.88359, 2128.89588, 2777.19525, 3783.07444])
        drag_v5_list_m = np.array([0, 8.47033, 31.64271, 68.48869, 118.70964, 182.09595, 258.82807, 349.47855, 456.92070, 550.00000, 650.00000, 918.59806, 1071.00961, 1245.39758, 1816.44728, 2369.96330, 3936.75503])
        drag_v5_list = [i*1 for i in drag_v5_list_m]
        v_list = np.array([0, 33, 66, 99, 132, 165, 198, 231, 264, 280, 310, 345, 360, 400, 500, 600, 800, 1200])

        drag_v6_list = [0, 17.463, 65.955, 143.612, 248.701, 382.638, 551.93, 664.104]
        v_v6_list = [0, 50, 100, 150, 200, 250, 300, 325, 343, 380, 430, 500, 600, 700, 800]

        cd_list = []
        for x in range(0,len(v_list)-1,1):
            cd = ( drag_v5_list[x] ) / ( 0.5 * rho_0 * pi * 3.25**2 * 0.0254**2 * ( v_list[x] + 0.001 )**2 )
            cd_list.append(cd)

        mach_list = v_list / solidworks_sos
        re_list = ( v_list * rocket_length * rho_0 ) / solidworks_mu
        current_mach = velocity / altitude_mach #mach at any given point, ranges from 0 - 2.3ish

        if current_mach < 0.3: #defined as the mach at which the mach number plays more role than re
            altitude_mu = mu_0*(altitude_temperature/t_0)**1.5*(t_0+s_0)/(altitude_temperature+s_0)
            current_re = velocity * altitude_density * rocket_length / altitude_mu
            for i in range(0,len(re_list),1):
                if current_re >= re_list[i] and current_re <= re_list[i+1]:
                    m = ( cd_list[i+1] - cd_list[i] ) / ( re_list[i+1] - re_list[i] )
                    current_cd = ( current_re - re_list[i] ) * m + cd_list[i]
                    break
                else:
                    continue
        else:
            for i in range(0,len(mach_list),1):
                if current_mach >= mach_list[-2]:
                    m = ( cd_list[-1] - cd_list[-2] ) / ( mach_list[-1] - mach_list[-2] )
                    current_cd = ( current_mach - mach_list[-1] ) * m + cd_list[-1]
                    break
                elif current_mach >= mach_list[i] and current_mach <= mach_list[i+1]:
                    m = ( cd_list[i+1] - cd_list[i] ) / ( mach_list[i+1] - mach_list[i] )
                    current_cd = ( current_mach - mach_list[i] ) * m + cd_list[i]
                    break
                else:
                    continue

        # G=1.45e6; #Effective Shear Modulus
        # Length1=24
        # Length2=20
        # Thickness=.22
        # Base=8.4
        # S=(Length1+Length2)/2*Base
        # AR=Base**2/S
        # e = 0.7
        #
        # final_cd = current_cd + current_cl**2 / (3.1415*AR*e)

        return altitude_temperature, altitude_pressure, altitude_density, current_mach, current_cd

    def force_calc():

        force_gravity = ( ( gravitational_constant * (total_rocket_mass) * (earth_mass) ) / (mojave_radius + height)**2 )
        force_drag = ( 0.5 * current_cd * altitude_density * velocity**2 * rocket_area )

        return force_gravity, force_drag

    def flutter_check(altitude_temperature,altitude_pressure,current_mach):

        G=1.45e6; #Effective Shear Modulus
        Length1=24
        Length2=16.75
        Thickness=.21
        Base=5
        S=(Length1+Length2)/2*Base
        AR=Base**2/S
        Lam=Length2/Length1

        temp_fahren = ( altitude_temperature - 273.15 ) * 1.8 + 32
        pres_lbft2 = 0.000145038 * altitude_pressure
        A=math.sqrt(1.4*1716.59*(temp_fahren+460))
        v=A*math.sqrt((G)/((1.337*AR**3*pres_lbft2*(Lam+1))/(2*(AR+2)*(Thickness/Length1)**3)))
        mach_flutter = v/A

        if current_mach >= mach_flutter:
            mach_flutter_list.append(mach_flutter)
            print(mach_flutter, current_mach)

        chance_flutter = ( current_mach - mach_flutter ) / mach_flutter

    def cg_v_t():

        # All of these calculations are in IPS, and should be kept that way
        # returns center of gravity in inches from top of nose cone

        if time >= burn_time:
            return dry_rocket_cg
        else:
            current_fuel_mass = fuel_mass_lb - (time/burn_time)*fuel_mass_lb
            fuel_mass_cg = (tank_bottom_position) - (tank_height*(1 - (time/burn_time)))/2
            current_cg = (current_fuel_mass*fuel_mass_cg + dry_rocket_mass_lb*dry_rocket_cg)/(current_fuel_mass+dry_rocket_mass_lb)
            return current_cg

    def rail_velocity():

        for i in range(0,len(height_list),1):
            deltax = rail_height - height_list[i]
            if deltax < 0:
                break

        return velocity_list[i]

    def sound_acceleration():

        for i in range(0,len(current_mach_list),1):
            if current_mach_list[i] > 0.8 and current_mach_list[i] < 1.2:
                if acceleration_list[i] > 0:
                    s_acc_list.append(acceleration_list[i])

    def center_of_pressure():

        # All of these calculations are in inches, and should be kept that way

        min_margin = 1.2 # calibers of stability

        # numbers in inches, for nose cone dimensions (VON KARMAN)
        nose_l = 32.5
        func1 = lambda x: (3.25*sqrt(math.acos(1-((2*x)/32.5))-(sin(2*math.acos(1-((2*x)/32.5)))/2)))/sqrt(pi) #integrate Von Karman, find area
        func2 = lambda x: x*(3.25*sqrt(math.acos(1-((2*x)/32.5))-(sin(2*math.acos(1-((2*x)/32.5)))/2)))/sqrt(pi) #Find centroid
        nose_a = (integrate.quad(func1, 0, 32.5))*2 #find nose cone total projected area
        pnose_d = (integrate.quad(func2, 0, 32.5)) #weird intermediate step
        nose_d = pnose_d[0]/nose_a[0] #nose cone centroid distance from nose cone tip (origin)

        # numbers in inches, for nose cone dimensions (POWER SERIES n=0.64)
        # nose_l = 32.5
        # func1 = lambda x: (3.25*(x/32.5)**0.64) #integrate Von Karman, find area
        # func2 = lambda x: x*(3.25*(x/32.5)**0.64) #Find centroid
        # nose_a = (integrate.quad(func1, 0, 32.5))*2 #find nose cone total projected area
        # pnose_d = (integrate.quad(func2, 0, 32.5)) #weird intermediate step
        # nose_d = pnose_d[0]/nose_a[0] #nose cone centroid distance from nose cone tip (origin)

        # numbers in inches, for rocket body dimensions
        body_diam = 6.5
        body_l = rocket_length-nose_l
        body_a = body_diam*(body_l) #projected body area
        body_d = nose_l+(body_l/2) #body centroid distance

        # numbers in inches, for fin dimensions
        rchord_l = 24 # dimensions for trapezoidal fin
        tchord_l = 16.75
        fheight = 4.64
        sweep_l = 6.612 # root chord segment for 55 degree angle of leading edge

        center_of_pressure = dry_rocket_cg + min_margin*rocket_radius_in*2

        fint_a = ((fheight*sweep_l)/2)*2 #first fin triangle area
        fint_d = nose_l+body_l-rchord_l+sweep_l/2 #centroid first fin triangle
        finr_a = (fheight*(rchord_l-sweep_l))*2 #rectangular fin portion area
        finr_d = nose_l+body_l-rchord_l+sweep_l+(rchord_l-sweep_l)/2 #centroid rectangular fin area
        # finbt_a = (fheight/2)*2 #last small triangle on fin area
        # finbt_d = nose_l+body_l-((rchord_l-tchord_l-sweep_l)*2)/3 #last small traingle centroid distance

        x = Symbol('x')
        var = solve((body_a+nose_a[0]+fint_a+finr_a+2*x*fheight)*center_of_pressure-(nose_l*nose_a[0]+body_a*body_d+fint_a*fint_d+finr_a*finr_d+2*x*fheight*(rocket_length+x)),x)

        # (nose_a[0]*nose_d+body_a*body_d+fint_a*fint_d+finr_a*finr_d+finbt_a*finbt_d)/(nose_a[0]+body_a+fint_a+finr_a+finbt_a)
        return var, center_of_pressure

    def append_lists():

        acceleration_list.append(acceleration_y)
        velocity_list.append(velocity)
        velocity_x_list.append(velocity_x)
        velocity_y_list.append(velocity_y)
        height_list.append(height)
        displacement_list.append(displacement)
        time_list.append(time)

        drag_coefficient_list.append(current_cd)
        current_mach_list.append(current_mach)
        force_gravity_list.append(force_gravity)
        force_drag_list.append(force_drag)
        density_list.append(altitude_density)
        temperature_list.append(altitude_temperature)
        pressure_list.append(altitude_pressure)
        cg_list.append(current_cg)

        theta_list.append(theta)
        omega_list.append(omega)
        alpha_list.append(alpha)

    def plot_plots():

        # theta_list_new = [i*(180/pi) for i in theta_list]
        # plt.plot(time_list,theta_list_new,'r')
        # plt.xlabel('Time (s)')
        # plt.ylabel('Angle of Attack (degrees)')
        # plt.show()

        color_list = []

        for x in range(0,6,1):
            color = "%06x" % random.randint(0, 0xFFFFFF)
            color_2 = '#' + color
            color_list.append(color_2)

        ft_list = np.asarray(height_list) * meters_to_feet
        plt.subplot(3,2,1)
        plt.plot(time_list, ft_list, color_list[0])
        plt.ylabel('Height (ft)')
        plt.suptitle('"Baby Come Back"',fontsize=16)

        fts_list = np.asarray(velocity_list) * meters_to_feet
        plt.subplot(3,2,2)
        plt.plot(time_list, fts_list, color_list[1])
        plt.ylabel('Velocity (ft/s)')

        ftss_list = np.asarray(acceleration_list) * meters_to_feet
        plt.subplot(3,2,3)
        plt.plot(time_list, ftss_list, color_list[2])
        plt.ylabel('Acceleration (ft/s^2)')
        plt.xlabel('Time (s)')

        # plt.subplot(3,2,4)
        # plt.plot(time_list, cg_list, color_list[3])
        # plt.ylabel('Center of Gravity (in from N.C.)')
        # plt.xlabel('Time (s)')

        plt.subplot(3,2,4)
        plt.plot(time_list, velocity_x_list, color_list[3])
        plt.ylabel('Velocity (x) m/s')
        plt.xlabel('Time (s)')

        theta_list_new = [i*(180/pi) for i in theta_list]
        plt.subplot(3,2,5)
        plt.plot(time_list,theta_list_new,'r')
        plt.xlabel('Time (s)')
        plt.ylabel('Angle of Attack (degrees)')

        plt.subplot(3,2,6)
        plt.plot(time_list, displacement_list, color_list[3])
        plt.ylabel('Displacement (m)')
        plt.xlabel('Time (s)')

        plt.subplots_adjust(left=0.2,wspace=0.4,hspace=0.2,top=0.9)
        plt.show()

        # cgg_list=[]
        # for i in range(0,len(cg_list),1):
        #     cgg = (-cg_list[i] + 8.35*12)/6.5
        #     cgg_list.append(cgg)
        # plt.plot(time_list, cgg_list, 'r')
        # plt.ylabel('Calibers of Stability')
        # plt.xlabel('Time (s)')
        # plt.show()

    """ Fundamental Constants """

    gravitational_constant = 6.67408e-11 #constant capital G
    avg_g = 9.81 #just used for calculation purposes
    boltzmann_constant = 1.38064852e-23 #in units of m^2 kg s^-2 K^-1
    gamma = 1.4
    gas_constant = 287.05 #J/(kg*K)
    gas_R = 8.314 #the other gas constant
    solidworks_sos = 345.45 #calculated with values, m/s
    solidworks_mu = 1.8213e-5 #SolidWorks viscosity at 297K
    mu_0 = 1.716e-5 #kg/m-s, standard values`
    t_0 = 273.11 #K, used for altitude_viscosity calc
    s_0 = 110.56 #K, used for altitude_viscosity calc
    rho_0 = 1.1826 #kg/m^3, from SolidWorks
    meters_to_feet = 3.28084 #meter to feet conversion
    lbf_to_n = 4.4482216152605 #pound force to newton conversion
    lb_to_kg = 0.453592
    pa_to_psi = 0.000145038
    real_e = 2.71 #approximate value for e
    sealevel_pressure = 101325 #atmospheric pressure at sea level in Pascals
    earth_mass = 5.972e24 #kilograms
    M = 0.0289644 #kg/mol, average mass of atmosphere
    temp_slope = 2 / 1000 #temperature decreases by 2 degrees for every increase in 1000ft altitude
    dt = 0.01 #increments of time value
    dt2 = 0.05

    """ Mojave Desert Specific Values """

    mojave_radius = 6371.13738e3 #distance in meters from mojave latitude (35 degrees) to center of earth
    mojave_temperature = 23.8889 + 273.15 #degrees kelvin, approximate for May 2020
    mojave_pressure = 100846.66 #pressure in Pascals, converted from mmHg
    mojave_density = mojave_pressure / ( gas_constant * mojave_temperature )
    avg_airmass = 28.95 / 6.022e20 #average mass of a single air molecule in kg
    max_wind_speed = 13.41 #wind speed in m/s, equivalent to 30mph
    rail_height = 60 / meters_to_feet #rail height in meters

    """ Rocket Constants """

    #converts all inputs from imperial to metric
    #thrust = thrust_lbf * lbf_to_n #convets thrust to N for calculations
    fuel_mass = fuel_mass_lb * lb_to_kg
    dry_rocket_mass = dry_rocket_mass_lb * lb_to_kg

    burn_time = 9208 / ( thrustofAlt(0) * 1/lbf_to_n) #estimated total burn time of liquid fuel, in seconds
    total_impulse = 9208 #lbf*s, specified by FAR MARS as total allowable impulse
    total_rocket_mass = fuel_mass + dry_rocket_mass #total mass, in kg
    mass_change = ( fuel_mass / burn_time ) * dt #assuming constant change of mass

    rocket_radius = ( innerRadius_in / 12 ) * ( 1 / meters_to_feet ) #radius of the rocket, meters
    rocket_area = pi * ( rocket_radius**2 ) #area from radius
    rocket_length = 11.92*12 #total length, from nose cone to bottom sleeve, in inches
    rocket_roughness = 3e-6 #surface roughness of carbon fiber in meters
    dry_rocket_cg = 7.07*12 #from feet to inches, measured from top of nose cone
    tank_bottom_position = 9.85*12 #distance from top of nose in inches
    tank_height = 4*12 #tank height in inches

    """ Initialize Variables """

    time = 0.0 #sets all variables to zero for start of the flight
    height = 0.0
    displacement = 0.0
    velocity_x = 0.0
    velocity_y = 0.0
    acceleration_x = 0.0
    acceleration_y = 0.0

    height_track = 0.0
    velocity_track = 0.0
    mach_track = 0.0
    q_track = 0.0
    flutter_track = 0.0
    current_impulse = 0.0

    theta = 0
    omega = 0
    alpha = 0
    current_torque = 0
    rotation_time = 0

    acceleration_list = [] #sets empty lists to be filled with data
    velocity_list = []
    velocity_x_list = []
    velocity_y_list = []
    height_list = []
    displacement_list = []
    time_list = []

    drag_coefficient_list = []
    current_mach_list = []
    force_gravity_list = []
    force_drag_list = []
    density_list = []
    temperature_list = []
    pressure_list = []
    mach_flutter_list = []
    s_acc_list = []
    cg_list = []

    theta_list = []
    omega_list = []
    alpha_list = []
    rotation_time_list = []

    while current_impulse <= total_impulse:

        if ((len(time_list)-1) % 50) == 0:
            wind_speed = max_wind_speed*random.random()

        # wind_speed = max_wind_speed

        velocity = sqrt(velocity_x**2+velocity_y**2)

        altitude_temperature, altitude_pressure, altitude_density, current_mach, current_cd = cd_altitude_calc() #finds density and drag coefficient at any given iteration
        force_gravity, force_drag = force_calc() #finds gravity and drag at any time
        flutter_check(altitude_temperature,altitude_pressure,current_mach)
        current_cg = cg_v_t()
        thrust = thrustofAlt(altitude_pressure) #1000 * thrust_change()

        if height <= rail_height:
            acceleration_y = ( thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
        else:
            theta, omega, alpha = rotation(wind_speed,theta)
            acceleration_y = ( thrust * cos(theta) - force_gravity - force_drag * cos(theta) ) / total_rocket_mass #finds y component acceleration
            acceleration_x = ( thrust * sin(theta) - force_drag * sin(theta) ) / total_rocket_mass #finds x component of acceleration

            # acceleration_x = ( thrust * sin(theta) - force_drag * sin(theta) )
            # if thrust < ( force_drag + force_gravity ):
            #     force_drag = thrust - force_gravity
            #     acceleration = ( thrust * cos(theta*(180/pi)) - force_gravity - force_drag * cos(theta*(180/pi)) ) / total_rocket_mass #finds resultant acceleration
            #     print("Uh oh. Rocket's going down")

        velocity_y += ( acceleration_y * dt ) #increment vertical velocity by small steps
        velocity_x += ( acceleration_x * dt ) #increment horizontal velocity
        # velocity_x = tan(theta) * velocity_y #increment horizontal velocity by small steps
        height += ( velocity_y * dt ) #same process for height
        displacement += ( velocity_x * dt ) #same process for displacement
        time += dt #increase total time by dt
        total_rocket_mass -= mass_change #at this time the rocket is losing fuel mass
        current_impulse += (thrust/lbf_to_n) * dt

        omega += alpha*dt
        theta -= omega*dt

        if velocity > velocity_track: #keeps track of the greatest total speed
            altitude_mach = sqrt( gamma * gas_constant * altitude_temperature ) #mach from temperature
            velocity_track = velocity
            mach_track = velocity / altitude_mach

        q = velocity**2 * altitude_density * 0.5
        if q > q_track:
            q_track = q

        append_lists() # calculation party time
    burn_time = time
    if(bTimeCalc):
        return height,time
    while velocity_y > 0:

        if ((len(time_list)-1) % 50) == 0:
            wind_speed = max_wind_speed*random.random()

        # wind_speed = max_wind_speed

        velocity = sqrt(velocity_x**2+velocity_y**2)

        altitude_temperature, altitude_pressure, altitude_density, current_mach, current_cd = cd_altitude_calc()
        force_gravity, force_drag = force_calc()
        flutter_check(altitude_temperature,altitude_pressure,current_mach)
        current_cg = cg_v_t()
        theta, omega, alpha = rotation(wind_speed,theta)

        acceleration_y = ( - force_gravity - force_drag * cos(theta) ) / total_rocket_mass
        acceleration_x = ( - force_drag * sin(theta) ) / total_rocket_mass
        velocity_y += ( acceleration_y * dt2 )
        velocity_x += ( acceleration_x * dt2 )
        # velocity_x = tan(theta) * velocity_y #increment horizontal velocity by small steps
        height += ( velocity_y * dt2 )
        displacement += ( velocity_x * dt2 ) #same process for displacement
        time += dt2

        omega += alpha*dt
        theta -= omega*dt

        if height > height_track: #stops the simulation when rocket is at apogee
            height_track = height

        append_lists() # calculation party time pt. 2

#    print('')
#    print("Max velocity =", round(velocity_track*meters_to_feet,3), "ft/s")
#    print("Max mach =", round(mach_track,3))
#    print("Max dynamic pressure =", round(q_track*pa_to_psi,3), "psi")
#    print("Max drag force =", round(max(force_drag_list)/lbf_to_n,3), "lbf")
#    print("Max acceleration =", round(max(acceleration_list)*meters_to_feet,3), "ft/s^2; (", round(max(acceleration_list)/9.8,3), "g's )")
#    var, cop = center_of_pressure()
#    print("Center of pressure =", round(cop/12,3), "ft from top, Extension = [", round(var[1],3), "inches ]")
#    print("Velocity off the rail =", round(rail_velocity()*meters_to_feet,3), "ft/s")
#    sound_acceleration()
#    print("Acceleration through transonic = [", round(s_acc_list[0]*meters_to_feet,3),",", round(s_acc_list[-1]*meters_to_feet,3),"] ft/s^2")
#    print("Center of Gravity = [ Worst:", round(max(cg_list)/12,3),", Initial:", round(cg_list[0]/12,3), "Final:", round(cg_list[-1]/12,3),"] ft from N.C.")
#    print("Caliber of Stability = [ Min:", round((-max(cg_list)+cop)/(6.5),3),", Max:", round((-min(cg_list)+cop)/(6.5),3),"]")
#    print("Maximum angle of attack = [", round(max(theta_list)*(180/pi),3), "degrees ]")
#    if len(mach_flutter_list) >= 1:
#        print("Fin Status: [Removed]")
#    else:
#        print("Fin Status: [Attached]")
#    print("Apogee =" , round(height_track*meters_to_feet,3), "ft")
#    print('')
    return round(height_track*meters_to_feet,3), burn_time
    #plot_plots()

# thrust_lbf = float(args[1])
# rocket_radius_in = float(args[2])
# fuel_mass_lb = float(args[3])
# dry_rocket_mass_lb = float(args[4])
#def thrust_change(height):
#    
#    alt_km = [0.0000, 0.2730, 0.5480, 0.8220, 1.0960, 1.3700, 1.6440, 1.9180, 2.1920, 2.4670, 2.7410, 3.0150, 3.2890, 3.5630, 3.8370, 4.1110, 4.3850, 4.6590, 4.9330, 5.2070, 5.4810]
#    thrust_kn = [3.1630, 3.1760, 3.1880, 3.2000, 3.2120, 3.2240, 3.2350, 3.2460, 3.2570, 3.2670, 3.2770, 3.2870, 3.2960, 3.3060, 3.3150, 3.3230, 3.3320, 3.3400, 3.3480, 3.3560, 3.3640]
#    sir_thrust = 0
#    for i in range(0,len(alt_km),1): #creates a linear regression for the above data, with individual equations for lines connecting each of the points
#        if height/1000 > alt_km[-1]: #if height is greater than 5.481km, assume linearly increasing thrust
#            m1 = (thrust_kn[-1] - thrust_kn[-2])/(alt_km[-1] - alt_km[-2])
#            sir_thrust = (height/1000-alt_km[-1]) * m1 + thrust_kn[-1]
#            print(sir_thrust)
#            break
#        elif height/1000 >= alt_km[i] and height/1000 < alt_km[i+1]:
#            m2 = (thrust_kn[i+1] - thrust_kn[i])/(alt_km[i+1] - alt_km[i])
#            sir_thrust = (height/1000-alt_km[i]) * m2 + thrust_kn[i]
#            break
#
#    return sir_thrust * 1000
 #determines initial parameters, can be changed
rocket_radius_in = 3.25
fuel_mass_lb = 35 #defined from the RPL spec sheet
dry_rocket_mass_lb = 59

# for z in range(0,100,1):
#print(main(thrust_change,rocket_radius_in,fuel_mass_lb,dry_rocket_mass_lb))
