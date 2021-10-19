import numpy as np
import pandas as pd
import pdb
import os
from modules.NoiseCraft import NoiseCraft
def main():
    
    equip = 'A320-211' # can be improved
    folder = 'AMS_LGW'
    #folder = 'ORY_to_BIA'
    meteo_fn = 'meteo.csv'
    flight_id = '520812451'
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
    
    #craft.calculate_distance(origin.to_numpy(), column_names)
    #craft.lat_lon_to_m(origin.to_numpy(), column_names)
    #craft.calculate_CAS(column_names)
    
    
    # TEMPORARY STEP SINCE HAVERSINE IS NOT PRECISE
    radar_aug = pd.read_csv(os.path.join('data', 'radar_tracks_n_meteo', 
        folder, 'flight_augmented.csv'))
    craft.Radar = radar_aug
    
    craft.load_data(column_names)
    
    mincount = 3
    maxcount = 5

    after_TO = True
    
    #Not working wrt to distance. Numpy error??
    craft.segment(mincount, maxcount, wrt_to_time=True, after_TO=after_TO, 
            normalize = True)
    
    craft.recognize_steps()
    craft.plot_segmented()
    if after_TO:
        craft.extrapolate_TO_distance()
    print('Steps: ' + str(craft.steps))
    craft.plot_segmented()

    craft.new_vert_profile(column_names)
    return 0, craft, origin


if __name__ == '__main__':
    success, craft, ori = main()
    print('Exited: ' + str(bool(~success)))
