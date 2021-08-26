import pandas as pd
import os
class NoiseCraft:
    
    def __init__(self, ACFT_ID):
        self.ACFT_ID = ACFT_ID
        path = r'./data/ANP_Database/'
        dirs = ['Aerodynamic_coefficients.csv', 'Aircraft.csv', 
                'Default_approach_procedural_steps.csv', 
                'Default_departure_procedural_steps.csv', 'Default_weights.csv',
                'Jet_engine_coefficients.csv']
        #Only added for A320-211, categories missing like propeller engine,etc.
        self.Aerodynamic_coefficients = pd.read_csv(os.path.join(path, dirs[0]),
                 sep=';')
        self.Aircraft = pd.read_csv(os.path.join(path, dirs[1]), sep=';')
        self.Approach_steps = pd.read_csv(os.path.join(path, dirs[2]), sep=';')
        self.Departure_steps = pd.read_csv(os.path.join(path, dirs[3]), sep=';')
        self.Default_weights = pd.read_csv(os.path.join(path, dirs[4]), sep=';')
        self.Jet_engine_coefficents = pd.read_csv(os.path.join(path, dirs[5]),
                sep=';')

