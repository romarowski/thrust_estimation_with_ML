import numpy as np
import pdb
def pressure_ratio_(Alt, Papt):
    # Papt [inHg], pressure at airport 
    # Alt [ft], aircraft altitude above airport
    
    g = 32.17 #ft/s^2 gravity
    B = 0.003566 #°F/ft atmospheric lapse
    R = 1716.59 #ft*lb/(sl*°R)
    T0 = 518.67 #°R temperature at std mean sea level
    P0 = 29.92 #inHg pressure at std mean sea level
    gBR = np.round(g/B/R, 4)
    
    np.seterr('raise')
    delta = ((Papt / P0)**(1/gBR) - (B * Alt / T0))**gBR

    return delta

def temperature_ratio_(Alt, Tapt):
    # Tapt [°F], temperature at airport
    # Alt [ft], aircraft altitude above airport
  
    B = 0.003566 #°F/ft atmospheric lapse
    FtoR = 459.67 # Conversion factor from °F to °R
    T0 = 518.67 #°R temperature at std mean sea level

    Theta = (FtoR + Tapt - B * Alt) / T0 

    return Theta

def air_density_ratio(Alt, Tapt, Papt):
    # Tapt [°F], temperature at airport
    # Papt [inHg], pressure at airport 
    # Alt [ft], aircraft altitude above airport
    delta = pressure_ratio_(Alt, Papt)
    Theta = temperature_ratio_(Alt, Tapt)

    sigma = delta / Theta

    return sigma