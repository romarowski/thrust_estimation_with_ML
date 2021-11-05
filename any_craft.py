import numpy as np
import pandas as pd
import pdb
import os
from modules.NoiseCraft import NoiseCraft
def main():
    
    plot = [] 
    while plot != 'y'  and plot != 'n':
        plot = input('plot? (y/n) ')
        plot_ = plot == 'y'

    
    
    equip = 'A320-211' # can be improved
    
    folder = 'AMS_LGW'
    flight_id = '520812451'
    
   #folder = 'ORY_to_BIA'
    
    #folder = 'AMS_TNG'
    #flight_id = '520831277'
    
    meteo_fn = 'meteo.csv'
    
    
    date = '20190101'
    radar_fn  = date + '_' + flight_id + '.csv'
    craft = NoiseCraft(equip)
    craft.load_meteo(folder, meteo_fn)
    craft.load_radar(folder, radar_fn)
    is_departure = True
    column_names = {'Pressure'    : 'pressure (hpa)',
                    'WindDir'     : 'wind direction (deg)',
                    'WindSpeed'   : 'wind speed (kts)',
                    'Heading'     : 'heading',
                    'Altitude'    : 'altitude',
                    'Temperature' : 'temperature (°C)',
                    'GroundSpeed' : 'speed',
                    'Lat'         : 'latitude',
                    'Lon'         : 'longitude'
                   }
    origin = craft.Radar.iloc[0][[column_names['Lat'], column_names['Lon']]]
    
    
    craft.lat_lon_to_m(origin.to_numpy(), column_names)
    craft.calculate_distance()
    craft.calculate_CAS(column_names)
    
    
    # TEMPORARY STEP SINCE HAVERSINE IS NOT PRECISE
    #radar_aug = pd.read_csv(os.path.join('data', 'radar_tracks_n_meteo', 
    #     folder, 'flight_augmented.csv'))
    #craft.Radar = radar_aug
    
    craft.load_data(column_names)
    
    mincount = 11 
    maxcount = 11

    after_TO = True
    
    craft.segment(mincount, maxcount, wrt_to_time=True, after_TO=after_TO, 
            normalize = False)
    
    craft.recognize_steps()


    
    if plot_:         
        craft.plot_segmented()
    
    if after_TO:
        craft.extrapolate_TO_distance()
    print('Steps: ' + str(craft.steps))
    
    #Look for point of thrust cutback
    #craft.segment(2, 2, wrt_to_time=True, thrust_cutback = True)
    
    
    #craft.time_seg = [  0.        ,  66.        ,  85.00230524,153, 183.1,
    #    8968075, 255.85383434, 303.69500753, 330.        ]
    #craft.cas_seg = [  0.        , 141.13015635, 143.74079402,249, 284.84,
    #    466217,      301.87449528, 303.77678871, 303.7783297 ]
    #craft.h_seg = [   0.        ,    0.        ,  887.50646384,2000, 3037,
    # .94248637, 6277.46698913, 8416.11229019, 9923.3773426 ]

    #craft.d_seg = [     0.        ,   7323.89550353,  12130.77228776, 344,
    #    80, 47415.8105425 , 84803.35358334, 112137.13826542, 126470.04953434]




    craft.new_vert_profile(column_names)
    #craft.new_vert_profile_pinv(column_names)
    
    if plot_:
        craft.plot_segmented()
        craft.plot_ANP_profile(column_names)

    craft.map_flaps()
    craft.generate_ANP_user_steps()
    craft.gen_fixed_point_profiles()

    return 0, craft, origin


if __name__ == '__main__':
    success, craft, ori = main()
    print('Exited: ' + str(bool(~success)))
