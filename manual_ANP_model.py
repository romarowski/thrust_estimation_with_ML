import os
import pandas as pd
from modules import NoiseCraft
from modules import ISA
import numpy as np
import pdb
from modules.ANP_steps import sin_of_climb_angle, TO_ground_roll_distance,\
        first_climb_CAS, distance_accelerate, corrected_net_thrust

def mid_values(vector):
    return (vector[1:] + vector[:-1]) / 2
def procedural_steps_takeoff(segmented_input, acraft, treshold_cas=10):
    dist = segmented_input[0]
    alt  = segmented_input[1]
    cas  = segmented_input[2]
    cas_derivative = np.diff(cas) <= treshold_cas
    steps = ['Climb' if der else 'Accelerate' for der in cas_derivative]
    climbs_deltaH = np.diff(alt)[cas_derivative]
    climbs_deltaS  = np.diff(dist)[cas_derivative]
    sins_gamma = climbs_deltaH / np.sqrt(climbs_deltaH**2 + climbs_deltaS**2)
    return sins_gamma

def main():
    craft = NoiseCraft.NoiseCraft('A320-211')
    stage_length = 2
    weight = craft.Default_weights[craft.Default_weights['Stage Length']==\
            stage_length]['Weight (lb)'].iloc[0]
    weight = 144921.86326854
    path_ = os.path.join('radar_tracks', 'ORY_to_BIA')
    craft.load_meteo(os.path.join(path_, 'meteo.csv'))
    column_names = {'Pressure'    : 'press (mbar)',
                    'WindDir'     : 'wind dir (from)',
                    'WindSpeed'   : 'wind speed (kts)',
                    'Heading'     : 'heading',
                    'Altitude'    : 'altitude (ft)',
                    'Temperature' : 'temp (°C)',
                    'GroundSpeed' : 'ground speed (kts)',
                    'Lat'         : 'latitude',
                    'Lon'         : 'longitude'
                   }
    
    Tapt = craft.Meteo[column_names['Temperature']].iloc[0]
    Papt = craft.Meteo[column_names['Pressure']].iloc[0]


    Aero_coeff = craft.Aerodynamic_coefficients
    dep_Aero_coeff = Aero_coeff[Aero_coeff["Op Type"] == "D"]
    Tapt = (Tapt * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)
    Papt = Papt / 33.864 #Conversion mbar to inHg (not ideal)
    
    dist = np.array([6746.26552269, 25498.04607828, 42533.60420978,  
        76004.85788356, 100928.13324184])
    alt = np.array([792.59920924,  3952.48861837,  4443.47893222, 
        9103.27665055, 10524.94704538])
    cas = np.array([150.32682737, 150.794697  , 229.65403539, 229.66665054,
       270.65863514])
    steps = ['Climb', 'Accelerate', 'Climb', 'Accelerate']
    
    sigmas = ISA.air_density_ratio(alt, Tapt, Papt)
    mean_sigmas = mid_values(sigmas)
    mean_sigmas_ = ISA.air_density_ratio(mid_values(alt), Tapt, Papt)

    deltas = ISA.pressure_ratio_(alt, Papt)
    mean_deltas = mid_values(deltas)
    print(deltas)
    W_delta = weight * np.reciprocal(mean_deltas)
    
    first_climb = True

    dep_steps = craft.Departure_steps[['Step Number', 'Step Type', 
        'Thrust Rating', 'Flap_ID']]
    climb_Flap_ID = dep_steps[dep_steps['Step Type']==\
            'Climb']['Flap_ID'].unique()
    TO_Flap_ID = climb_Flap_ID[0] 
    N = craft.Aircraft['Number Of Engines'].iloc[0]

    Jet_eng_coeff = craft.Jet_engine_coefficients
    thrust_rating = {'Reg'   : {'C': 'MaxClimb', 'TO' : 'MaxTakeoff'},
                     'HiTemp': {'C': 'MaxClimbHiTemp', 'TO' : 'MaxTkoffHiTemp'}
                     }
    thrust_type = 'Reg'
    midpoint_Temps = ISA.temperature_profile(mid_values(alt), Tapt)
    sins_gamma = np.zeros(np.shape(steps))
    for i, step in enumerate(steps):
        if step == 'Climb' and first_climb:
            R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                    TO_Flap_ID]['R'].iloc[0]
               
            thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                    thrust_rating[thrust_type]['TO']].reset_index(drop=True)
            
            Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                    midpoint_Temps[i], mid_values(alt)[i])
            print(Fn_delta)
            sins_gamma[i] = sin_of_climb_angle(N, Fn_delta, W_delta[i], R,
                    mid_values(cas)[i])
            
            first_climb = False
        elif step == 'Climb' and first_climb==False:

            R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                    climb_Flap_ID[1]]['R'].iloc[0]
            thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                    thrust_rating[thrust_type]['C']].reset_index(drop=True)

            Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                    midpoint_Temps[i], mid_values(alt)[i])
   
            sins_gamma[i] = sin_of_climb_angle(N, Fn_delta, W_delta[i], R,
                    mid_values(cas)[i])
  
            
    segmented_input = (dist, alt, cas)
    sins_gamma_est = procedural_steps_takeoff(segmented_input, craft)
    print(sins_gamma)
    print(sins_gamma_est)




    return 0

if __name__ == '__main__':
    success = main()
    print(np.bool_(~success))
