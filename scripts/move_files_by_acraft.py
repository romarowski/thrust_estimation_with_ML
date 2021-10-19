import pandas as pd 
import os 
import numpy as np

type_ = 'A320'


data = pd.read_csv('20190101_flights.csv')
data = data[data['equip'] == type_].reset_index(drop=True)
data.to_csv('20190101_flights_'+type_, index=False)

list_fligghts = os.listdir('20190101_positions/')
aa = pd.read_csv('20190101_flights_+'type_+'A320.csv')
flight_ids = aa['flight_id'].to_numpy()


matching = np.array([fn for fn in list_fligghts if \
        any(flight in fn for flight in flight_ids)])


