import pandas as pd
import os
import numpy as np
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
    pass

    def load_meteo(self, filename):
        path = r'./data/'
        self.Meteo = pd.read_csv(os.path.join(path, filename))
    pass

    def load_radar(self, filename):
        path = r'./data/'
        self.Radar = pd.read_csv(os.path.join(path, filename))
    pass

    def calculate_CAS(self, column_names):
        if self.Meteo   = None:
            print('Load meteorological data to calculate CAS')
        elif self.Radar = None:
            print('Load radar/ADSB data to calculate CAS')
        else:
            #Solution for one windspeed/direction
            winddir    = self.Meteo[column_names['WindDir']].iloc[0]
            windspeed = self.Meteo[column_names['WindSpeed']].iloc[0]
            gamma = np.pi - np.deg2rad(winddir) +\
                    np.deg2rad(self.Radar[column_names['Heading']])
            TAS=(self.Radar[column_names['GroundSpeed']]-windspeed*np.cos(gamma))

self.Radar['CAS (kts)']
