import pandas as pd
import os
import pdb
import numpy as np
from modules.ISA import air_density_ratio, pressure_ratio_, pressure_altitude,\
        temperature_profile, temperature_ratio_
from modules import inm_map_projection
from scipy import interpolate
from modules.haversine_distance import haversine
import matplotlib.pyplot as plt
from scipy import optimize, interpolate
from modules.ANP_steps import sin_of_climb_angle, TO_ground_roll_distance,\
        first_climb_CAS, distance_accelerate, corrected_net_thrust



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
        self.Jet_engine_coefficients = pd.read_csv(os.path.join(path, dirs[5]),
                sep=';')
        self._Meto_loaded = False
        self.Meteo = None
        self_Radar_loaded = False
        self.Radar = None
        pass

    def load_meteo(self, folder, filename):
        path = os.path.join('data', 'radar_tracks_n_meteo', folder) 
        #Improve,  makes no sense like this, '/' is Linux only
        self.Meteo = pd.read_csv(os.path.join(path, filename))
        self._Meteo_loaded = True
        pass

    def load_radar(self, folder, filename):
        path = os.path.join('data', 'radar_tracks_n_meteo', folder)
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
            self.Radar['TAS (kts)'] = TAS

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

    #def calculate_distance():
            #TODO
    
   # def calculate_distance(self, origin, column_names):
   #     self.Radar['Distance (ft)'] = self.Radar[[column_names['Lat'],\
   #             column_names['Lon']]].apply(haversine, axis=1, 
   #             origin = origin)
   #     pass
        
    #def filter_data(runway_geo):
        #TODO
    
    
    def load_data(self, column_names):
        flight = self.Radar
        d = flight['distance (ft)'].to_numpy()
        h = flight[column_names['Altitude']].to_numpy()
        cas = flight['CAS (kts)'].to_numpy()
        time = flight['snapshot_id'].to_numpy()
        time -= time[0]
        #cas_o = cas[ydata>0]
        #xdata_o = xdata[ydata>0]
        #ydata_o = ydata[ydata>0]
        #plt.subplot(221)
        #plt.plot(d, h, '.')
        #plt.xlabel('distance (ft)')
        #plt.ylabel('altitude (ft)')
        #plt.subplot(223)
        #plt.plot(d, cas, '.')
        #plt.xlabel('distance (ft)')
        #plt.ylabel('cas (kts)')
        #plt.subplot(222)
        #plt.plot(time, h, '.')
        #plt.xlabel('time (s)')
        #plt.ylabel('altitude (ft)')
        #plt.subplot(224)
        #plt.plot(time, cas, '.')
        #plt.xlabel('time (s)')
        #plt.ylabel('cas (kts)')
        #plt.show()
        ## Save info from radar prolly should be broken down
        self.h   = np.float64(h)
        self.cas = cas
        self.d   = d
        self.time = time
        pass

    def segment(self, mincount, maxcount, wrt_to_time=False, after_TO=False,
            normalize=False, cas_weight = 1000):
        

        if wrt_to_time:
            X = np.copy(self.time)
        else:
            X = np.copy(self.d)
        
        Y = np.copy(self.h)
        cas = np.copy(self.cas)
        
        if normalize:
            norm_x = X.max()
            norm_y = Y.max()
            norm_cas = cas.max()

            X = np.divide(X, norm_x)
            Y =  np.divide(Y, norm_y)
            cas = np.divide(cas, norm_cas)
            
            cas_weight = 1.0

        if after_TO:
            X = X[Y>0]
            cas = cas[Y>0]
            Y = Y[Y>0]

        
        xmin = X.min()
        xmax = X.max()
        n = len(X)
        
        AIC = []
        r   = []
        regions = []
        i = 0

        for count in range(mincount, maxcount + 1): 
            seg = np.full(count - 1, (xmax - xmin) / count)

            px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
            mask = [[np.abs(X - x) < (xmax - xmin) * 0.1] for x in px_init]


            py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.1].mean()\
                    for x in px_init])
            pcas_init = np.array([cas[np.abs(X - x) < (xmax - xmin) * 0.1].mean()\
                    for x in px_init])


            def func(p, count):
                seg = p[:count - 1]
                py = p[count - 1:2*count]
                pcas = p[2*count:]
                px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
                return px, py, pcas

            def err(p):
                px, py, pcas = func(p, count)
                Y2 = np.interp(X, px, py)
                CAS2 = np.interp(X, px, pcas)
                if after_TO:
                    
                    penalty_FC = max(0, np.diff(pcas)[0])**2
                    #if np.diff(pcas)[0] > 0:
                    #    penalty_FC = np.diff(pcas)[0]*1e2
                    #else:
                    #    penalty_FC = 0
                    penalty_TO = 0
                else:
                    
                    penalty_FC = max(0, np.diff(pcas)[1])**2
                    #if np.diff(pcas)[1] > 0:
                    #    penalty_FC = np.diff(pcas)[1]*1e2
                    #else: 
                    #    penalty_FC = 0
                    penalty_TO = max(0, np.diff(py)[0])**2
                    #if np.diff(py)[0] > 0:
                    #    penalty_TO = np.diff(py)[0]*1e2
                    #else:
                    #    penalty_TO = 0
                if not np.all(np.diff(pcas)>=0):
                    penalty_CAS = 10e3
                else:
                    penalty_CAS = 0
                if not np.all(np.diff(py) >= 0):
                    penalty_y = 10e3
                else: 
                    penalty_y = 0
                return np.mean((Y - Y2)**2) +\
                       cas_weight*np.mean((cas - CAS2)**2) + penalty_CAS + \
                       penalty_FC + penalty_TO + penalty_y
        
            x0=np.r_[seg, py_init, pcas_init]
            r.append(optimize.minimize(err, x0=x0, method='Nelder-Mead',
                    options={'adaptive': False, 'fatol':1e-5}))
            AIC.append(n * np.log10(err(r[i].x)) + 4 * count)
            #BIC = n * np.log10(err(r.x)) + 2 * count * np.log(n)

            regions.append(count)

            i += 1
        AIC_min = min(AIC)
        min_index = AIC.index(AIC_min)
        r_min = r[min_index]
        reg_min = regions[min_index]
        print('AICs:' + str(AIC))
        print('Regions:' + str(regions))
        
        if wrt_to_time:
            ptimes, py, pcas = func(r_min.x, reg_min)
            if normalize:
                ptimes = np.multiply(ptimes, norm_x)
            px = np.interp(ptimes, self.time, self.d)
        else:        
            px, py, pcas = func(r_min.x, reg_min)
            if normalize:
                px = np.multiply(px, norm_x)
            ptimes = np.interp(px, self.d, self.time)

        if normalize:
            py = np.multiply(py, norm_y)
            pcas = np.multiply(pcas, norm_cas)
        
        self.d_seg = px
        self.h_seg = py
        self.cas_seg = pcas
        self.time_seg = ptimes
        self.optimizer_result = r_min
        print('Optimization finished')
        pass

    def plot_segmented(self):
        px = self.d_seg
        py = self.h_seg
        pcas = self.cas_seg
        ptimes = self.time_seg
        X = self.d
        Y = self.h
        CAS = self.cas
        times = self.time
        
        plt.subplot(221)
        plt.plot(X, Y, '.', label='Radar')
        plt.plot(px, py, '-or', label='Model')
        plt.xlabel('Dist [ft]')
        plt.ylabel('Alt [ft]')
        plt.legend()
        
        plt.subplot(223)
        plt.plot(X, CAS, '.', label='Radar')
        plt.plot(px, pcas, '-or', label='Model')
        plt.xlabel('Dist [ft]')
        plt.ylabel('CAS [kts]')
        plt.legend()

        plt.subplot(222)
        plt.plot(times, Y, '.', label='Radar')
        plt.plot(ptimes, py, '-or', label='Model')
        plt.xlabel('Time [s]')
        plt.ylabel('Alt [ft]')
        plt.legend()

        plt.subplot(224)
        plt.plot(times, CAS, '.', label='Radar')
        plt.plot(ptimes, pcas, '-or', label='Model')
        plt.xlabel('Time [s]')
        plt.ylabel('CAS [kts]')
        plt.legend()
        
        plt.show()

    def recognize_steps(self, treshold_CAS=10, after_TO=True):
        cas = self.cas_seg
        cas_derivative = np.diff(self.cas_seg) <= treshold_CAS
        steps = ['Climb' if der else 'Accelerate' for der in cas_derivative]
        self.steps = steps
        pass

    def extrapolate_TO_distance(self):
        dist = self.d_seg
        alt  = self.h_seg
        cas  = self.cas_seg
        steps = self.steps
        time = self.time_seg
        f = interpolate.interp1d(alt[:2], dist[:2], fill_value = 'extrapolate')
        steps_ = np.copy(steps)
        steps_ = np.r_[['TakeOff'], steps_]
        dist_ = np.copy(dist)
        alt_  = np.copy(alt)
        cas_  = np.copy(cas)
        time_  = np.copy(time)
        dist_[0] = f(0)
        alt_[0] = 0.
        time_ = np.r_[0., time_]
        dist_ = np.r_[0., dist_]
        alt_ = np.r_[0., alt_]
        cas_ = np.r_[0, cas_]

        self.d_seg   = dist_
        self.h_seg   = alt_
        self.cas_seg = cas_
        self.steps   = steps_ 
        self.time_seg    = time_
        
        pass
    
    def new_vert_profile(self, column_names, op_type = 'D'):
        
        stage_length = 2
        weight = self.Default_weights[self.Default_weights['Stage Length']==\
        stage_length]['Weight (lb)'].iloc[0]

        dist   = self.d_seg
        alt    = self.h_seg
        cas    = self.cas_seg
        times  = self.time_seg
        steps = self.steps
        print('Laying down new profile...')
        def mid_values(vec):
            return (vec[1:] + vec[:-1]) / 2
        def sins_gamma_estimation(dist, alt, cas):
            climbs_deltaH = np.diff(alt)
            climbs_deltaS = np.diff(dist)
            sins_gamma = climbs_deltaH /\
                    np.sqrt(climbs_deltaH**2 + climbs_deltaS**2)
            return sins_gamma
        Tapt = self.Meteo[column_names['Temperature']].iloc[0]
        Papt = self.Meteo[column_names['Pressure']].iloc[0]

        Aero_coeff = self.Aerodynamic_coefficients
        dep_Aero_coeff = Aero_coeff[Aero_coeff["Op Type"] == op_type]
        Tapt = (Tapt * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)
        Papt = Papt / 33.864 #Conversion mbar to inHg (not ideal)

        sigmas = air_density_ratio(alt, Tapt, Papt)
        tas = np.multiply(cas, np.reciprocal(np.sqrt(sigmas)))
        mean_sigmas = mid_values(sigmas)
        mean_sigmas_ = air_density_ratio(mid_values(alt), Tapt, Papt)

        tas_geom_mean = np.sqrt(mid_values(np.power(tas, 2))) #Impact eq-17
        tas_diff = np.diff(np.power(tas, 2))
        deltas = pressure_ratio_(alt, Papt)
        mean_deltas = mid_values(deltas)
        #print(deltas)
        W_delta = weight * np.reciprocal(mean_deltas)

        
        
        first_climb = True


        #Obtain flaps by type of step: Accel, TO, Climb for now
        dep_steps = self.Departure_steps[['Step Number', 'Step Type',
            'Thrust Rating', 'Flap_ID']]
        climb_Flap_ID = dep_steps[dep_steps['Step Type']==\
                'Climb']['Flap_ID'].unique()
        TO_Flap_ID = dep_steps[dep_steps['Step Type']==\
                'Takeoff']['Flap_ID'].unique() #climb_Flap_ID[0] 
        accel_Flap_ID = dep_steps[dep_steps['Step Type']==\
                'Accelerate']['Flap_ID'].unique()

        N = self.Aircraft['Number Of Engines'].iloc[0]

        Jet_eng_coeff = self.Jet_engine_coefficients
        thrust_rating = {'Reg'   : {'C': 'MaxClimb', 'TO' : 'MaxTakeoff'},
                         'HiTemp': {'C': 'MaxClimbHiTemp', 'TO' : 'MaxTkoffHiTemp'}
                         }
        thrust_type = 'Reg'
        midpoint_Temps = temperature_profile(mid_values(alt), Tapt)

        sins_gamma = np.zeros(np.shape(steps))
        segment_accel = np.zeros(np.shape(steps))
        ROCs = np.multiply(np.diff(alt), np.reciprocal(np.diff(times)))

        for i, step in enumerate(steps):
            if step == 'TakeOff':
                
                B8 = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        TO_Flap_ID[0]]['B'].iloc[0]
                theta = temperature_ratio_(Alt=0, Tapt=Tapt)
                
                W_delta_TO = weight / pressure_ratio_(Alt=0, Papt=Papt)
                
                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['TO']].reset_index(drop=True)
                
                #Should be CAS of first point
                Fn_delta = corrected_net_thrust(thrust_coeff, cas[1],
                    temperature_profile(Alt=0, Tapt=Tapt), Alt=0, Papt=Papt)
                est_TO8 = TO_ground_roll_distance(B8, theta, W_delta_TO, N, 
                        Fn_delta)
            
            elif step == 'Climb' and first_climb:
                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        TO_Flap_ID[0]]['R'].iloc[0]

                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['TO']].reset_index(drop=True)

                Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                #print(Fn_delta)
                sins_gamma[i] = sin_of_climb_angle(N, Fn_delta, W_delta[i], R,
                        mid_values(cas)[i])

                #FIRST CAS

                C = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        TO_Flap_ID[0]]['C'].iloc[0]

                est_first_climb_CAS = first_climb_CAS(C, weight)
                radar_first_CAS = cas[i]


                first_climb = False


            elif step == 'Climb' and first_climb==False:

                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        climb_Flap_ID[1]]['R'].iloc[0]
                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['C']].reset_index(drop=True)

                Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)

                sins_gamma[i] = sin_of_climb_angle(N, Fn_delta, W_delta[i], R,
                        mid_values(cas)[i])
            elif step == 'Accelerate':

                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        accel_Flap_ID[1]]['R'].iloc[0]

                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                    thrust_rating[thrust_type]['C']].reset_index(drop=True)

                Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                    midpoint_Temps[i], mid_values(alt)[i], Papt)
                segment_accel[i] = distance_accelerate(tas_diff[i],
                    tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
        
        sins_gamma_est = sins_gamma_estimation(dist, alt, cas)

        print(steps)
        print('ANP climb angle: '  + str(sins_gamma))
        print('RADAR climb angle ' +  str(sins_gamma_est))
        print('ANP accel segment '+ str(segment_accel))
        print('RADAR accel segment ' + str(np.diff(dist)))
        print('ANP first CAS ' + str(est_first_climb_CAS))
        print('RADAR first CAS ' + str(radar_first_CAS))
        print('ANP TO8 ' + str(est_TO8))
        print('RADAR TO8 ' + str(np.diff(dist)[0]))

    #def segment_time(self, dist, dist_full, times):
    #    ptimes = np.interp(dist, dist_full, times)
    #    return ptimes 


