def He_mass(LOX_vol, CH4_vol, tank_pressure, He_pressure=6000):

    R=2077 #J kg/K, specific gas constant for He
    Z=1.18 #gHe compressibility at 6000 psi and 300K, from https://cds.cern.ch/record/1444601/files/978-1-4419-9979-5_BookBackMatter.pdf
    psi_to_Pa = 6894.76

    He_pressure *= psi_to_Pa
    initial_temperature = 300 #K
    He_pressure_losses = 50 #psi, arbitrary
    final_he_pressure = (tank_pressure+He_pressure_losses)*psi_to_Pa
    tank_pressure *= psi_to_Pa
    propellant_volume = LOX_vol+CH4_vol
    sigma = 2 * 10**9 #Pa, or 2 GPa, yeild stress estimate of carbon fiber
    rho = 1800 #kg/m^3, of carbon fiber
    
    #from pg 24 of "Cryogenic propellant tank pressurization" by Rob Hermsen
    
    m_He=((tank_pressure*propellant_volume)/(R*initial_temperature))*(1.66/(1-(Z*final_he_pressure/He_pressure)))#kg
    K = 2
    m_He*=K #kg

    #m0 is the starting mass of helium in the pressure vessel
    
    n_He=m_He*1000/4 #mols

    V_He = n_He*8.3145*initial_temperature/He_pressure #m^3

    #done with helium calculations, starting Tank mass calculations

    r = 0.01*2.54*3 #3 inches in m, this is the ID of the He tank
    V_cylinder = V_He - (4/3)*3.14159*r**3
    h = V_cylinder/(3.14159*r**2)
        
  
    wall_thickness = 1.5 * He_pressure * r / sigma #m
    #the spherical and cylindrical sections of a pressure vessel have diffent stresses, hence different thicknesses.
    #For reference, check out http://web.mit.edu/course/3/3.11/www/modules/pv.pdf
    #however, for our purposes we'll stay on the safe side and take the thicker one
    
    cylinder_vol = h * 3.14159 * ((r+wall_thickness)**2 - r**2) #m^3
    spherical_vol =(4/3)*3.14159*((r+wall_thickness)**3 - r**3)

    m_tank = rho * (cylinder_vol + spherical_vol)
    m_total = m_He + m_tank

    return m_total
        
