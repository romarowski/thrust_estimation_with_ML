from modules import NoiseCraft
import numpy as np
import os
def main():
    path_ = os.path.join('radar_tracks', 'ORY_to_BIA')
    flight = NoiseCraft.NoiseCraft('A320-211')
    flight.load_meteo(os.path.join(path_, 'meteo.csv'))
    flight.load_radar(os.path.join(path_, 'flight.csv'))
    # CAS_example
    #column_names = {'Pressure'    : 'Airport pressure (In-Hg)',
    #                'WindDir'     : 'Winddirection',
    #                'WindSpeed'   : 'Headwind (kts)',
    #                'Heading'     : 'Heading',   
    #                'Altitude'    : 'Altitude (ft)',
    #                'Temperature' : 'Airport temperature (°F)',
    #                'GroundSpeed' : 'Ground speed (kts)'
    #               }
    # ORY_to_BYA
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
    origin = np.array([48.723333, 2.379444]) #ORY
    flight.calculate_CAS(column_names)
    flight.lat_lon_to_m(origin, column_names)

    return flight


if __name__ == '__main__':
    flight = main()
