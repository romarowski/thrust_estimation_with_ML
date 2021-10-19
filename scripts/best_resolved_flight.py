import os
import pandas as pd
import numpy as np
import pdb

date = '20190101_'
dirx = os.path.join('data', 'radar_tracks', 'FR24') 
db_f = pd.read_csv(os.path.join(dirx,'20190101_flights_A320.csv')) 
allflights = db_f['flight_id'].to_numpy()
dist = float('inf')
file_ = ''
for flight in allflights:
    db = pd.read_csv(os.path.join(dirx, '20190101_positions_A320',
       date + str(flight) + '.csv'))
    dist_ = np.min(db[(db['altitude']<=1000) & (db['altitude'].diff()>0) &\
            db['altitude']>0]['altitude'])
    if dist_ < dist:
        dist = dist_
        file_ = flight
        db_min = db


