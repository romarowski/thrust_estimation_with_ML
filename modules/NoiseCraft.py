import pandas as pd
import os
import pdb
import numpy as np
from modules.ISA import air_density_ratio
from modules import inm_map_projection
class NoiseCraft:
    
    def __init__(self, ACFT_ID):
        self.ACFT_ID = ACFT_ID
        path = r'./data/ANP_Database/'
        dirs = ['Aerodynamic_coefficients.csv', 'Aircraft.csv', 
                'Default_approach_procedural_steps.csv', 
                'Default_departure_procedural_steps.csv', 'Default_weights.csv',
                'Jet_engine_coefficients.csv']
        #Working e.g. for A320-211, categories missing like propeller engine,etc.
        self.Aerodynamic_coefficients = pd.read_csv(os.path.join(path, dirs[0]),
                 sep=';')
        self.Aircraft = pd.read_csv(os.path.join(path, dirs[1]), sep=';')
        self.Approach_steps = pd.read_csv(os.path.join(path, dirs[2]), sep=';')
        self.Departure_steps = pd.read_csv(os.path.join(path, dirs[3]), sep=';')
        self.Default_weights = pd.read_csv(os.path.join(path, dirs[4]), sep=';')
        self.Jet_engine_coefficents = pd.read_csv(os.path.join(path, dirs[5]),
                sep=';')
        self._Meto_loaded = False
        self.Meteo = None
        self_Radar_loaded = False
        self.Radar = None
    pass

    def load_meteo(self, filename):
        path = r'./data/' 
        #Improve,  makes no sense like this, '/' is Linux only
        self.Meteo = pd.read_csv(os.path.join(path, filename))
        self._Meteo_loaded = True
    pass

    def load_radar(self, filename):
        path = r'./data/'
        self.Radar = pd.read_csv(os.path.join(path, filename))
        self._Radar_loaded = True
    pass

    def calculate_CAS(self, column_names):
        if self._Meteo_loaded   == False:
            print('Load meteorological data to calculate CAS')
        elif self._Radar_loaded == False:
            print('Load radar/ADSB data to calculate CAS')
        else:
            #Solution for one windspeed/direction
            #pdb.set_trace()
            winddir    = self.Meteo[column_names['WindDir']].iloc[0]
            windspeed = self.Meteo[column_names['WindSpeed']].iloc[0]
            Tapt = self.Meteo[column_names['Temperature']].iloc[0]
            Tapt = (Tapt * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)
            Papt = self.Meteo[column_names['Pressure']].iloc[0]
            Papt = Papt / 33.864 #Conversion mbar to inHg (not ideal)
            gamma = np.pi - np.deg2rad(winddir) +\
                    np.deg2rad(self.Radar[column_names['Heading']])
            TAS=(self.Radar[column_names['GroundSpeed']]-windspeed*np.cos(gamma))

            sigma = air_density_ratio(self.Radar[column_names['Altitude']], 
                                      Tapt, Papt) 
            self.Radar['CAS (kts)'] = np.multiply(TAS, np.sqrt(sigma))

    pass

    def lat_lon_to_m(self, origin, column_names, proj_type='inm'):
        if self._Radar_loaded == False:
            print('Load radar/ADSB data to calculate CAS')
        elif proj_type == 'inm':
            parameters = inm_map_projection.parameters(origin)
            project = inm_map_projection.project
            self.Radar['x (m)'] , self.Radar['y (m)'] =\
                    zip(*self.Radar[[column_names['Lat'], column_names['Lon']]]\
                    .apply(project, axis = 1, parameters = parameters))

    pass
    
    #def calculate_distance(self):
        #TODO
