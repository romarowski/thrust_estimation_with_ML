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
        first_climb_CAS, distance_accelerate, corrected_net_thrust,\
        sin_of_climb_angle_pinv, first_climb_CAS_pinv, distance_accel_pinv,\
        corrected_net_thrust_pinv




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

    def calculate_distance(self):
        mt2ft = 3.281
        data = self.Radar[['x (m)', 'y (m)']].copy()
        data['distance (ft)'] = np.zeros((len(data), 1))

        for row in np.arange(1, len(data)):
            dist_x = (data['x (m)'].iloc[row-1] - data['x (m)'].iloc[row])**2
            dist_y = (data['y (m)'].iloc[row-1] - data['y (m)'].iloc[row])**2
            dist = np.sqrt(dist_x + dist_y) * mt2ft
            data['distance (ft)'].iloc[row] = data['distance (ft)'].iloc[row-1] +\
                    dist

        self.Radar['distance (ft)'] = data['distance (ft)']

        pass



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
    
    
    #def thrust_cutback(self, wrt_to_time
    
    
    
    def segment(self, mincount, maxcount, wrt_to_time=False, after_TO=False,
            normalize=False, cas_weight = 1000, first_climb_weight=2, 
            thrust_cutback = False):
        

        if wrt_to_time:
            X = np.copy(self.time)
        else:
            X = np.copy(self.d)
        
        Y = np.copy(self.h)
        cas = np.copy(self.cas)
        
        first_accel = True
        if thrust_cutback:
            for i, step in enumerate(self.steps):
                if step == 'Accelerate' and first_accel and wrt_to_time:
                    x_start = self.time_seg[i]
                    x_end   = self.time_seg[i+2]
                    Y = Y[(X>=x_start) & (X<=x_end)]
                    cas = cas[(X>=x_start) & (X<=x_end)]
                    X = X[(X>=x_start) & (X<=x_end)]
                    first_accel = False
                    break
                elif step == 'Accelerate' and first_accel and not wrt_to_time:
                    x_start = self.d_seg[i]
                    x_end   = self.d_seg[i+2]
                    X = X[(X>=x_start) & (X<=x_end)]
                    Y = Y[(X>=x_start) & (X<=x_end)]
                    cas = cas[(X>=x_start) & (X<=x_end)]
                    first_accel = False
                    break
            mincount = 2
            maxcount = 2
                     

        if normalize:
            norm_x = X.max()
            norm_y = Y.max()
            norm_cas = cas.max()

            X = np.divide(X, norm_x)
            Y =  np.divide(Y, norm_y)
            cas = np.divide(cas, norm_cas)
            
            cas_weight = 1.0

        else:
            norm_x = 1.0
            norm_y = 1.0
            norm_cas = 1.0

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
                    penalty_TO = 0
                else:
                    penalty_FC = max(0, np.diff(pcas)[1])**2
                    penalty_TO = max(0, np.diff(py)[0])**2
                
                penalty_CAS = max(0, np.mean(np.diff(pcas)))**2
                penalty_y   = max(0, np.mean(np.diff(pcas)))**2 
                #if not np.all(np.diff(pcas)>=0):
                #    penalty_CAS = 1e10
                #else:
                #    penalty_CAS = 0
                #if not np.all(np.diff(py) >= 0):
                #    penalty_y = 10e3
                #else: 
                #    penalty_y = 0
                
                #First climb precision
               
                if thrust_cutback:
                    cost = np.mean((Y-Y2)**2)+cas_weight*np.mean((cas-CAS2)**2)+\
                        penalty_CAS + penalty_y
                else:
                    cost = np.mean((Y - Y2)**2) + \
                           first_climb_weight*\
                           np.mean((Y[Y<1000./norm_y] - Y2[Y<1000./norm_y])**2) + \
                           cas_weight*np.mean((cas - CAS2)**2) + penalty_CAS + \
                                   penalty_FC + penalty_TO + penalty_y
                return cost

            x0=np.r_[seg, py_init, pcas_init]
            r.append(optimize.minimize(err, x0=x0, method='Nelder-Mead',
                    #options={'adaptive': False, 'fatol':1e-5},
                    )
                 )
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
        print('Region: ' + str(reg_min))
        
        if thrust_cutback:
            pdb.set_trace()

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
        print('Segmentation finished')
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
        
        print('Laying down new profile...')
        costs = []
        def loss(x):
            weight = x[0]
            Tflex  = x[1]

            #weight = 127496.7607

            dist   = self.d_seg
            alt    = self.h_seg
            cas    = self.cas_seg
            times  = self.time_seg
            steps = self.steps
            #print('Laying down new profile...')
            
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

            Tflex = (Tflex * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)

            
            Aero_coeff = self.Aerodynamic_coefficients
            dep_Aero_coeff = Aero_coeff[Aero_coeff["Op Type"] == op_type]
            Tapt = (Tapt * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)
            Papt = Papt / 33.864 #Conversion mbar to inHg (not ideal)
            

            sigmas       = air_density_ratio(alt, Tapt, Papt)
            tas          = np.multiply(cas, np.reciprocal(np.sqrt(sigmas)))
            mean_sigmas  = mid_values(sigmas)
            mean_sigmas_ = air_density_ratio(mid_values(alt), Tapt, Papt)

            tas_geom_mean = np.sqrt(mid_values(np.power(tas, 2))) #Impact eq-17
            tas_diff      = np.diff(np.power(tas, 2))
            deltas        = pressure_ratio_(alt, Papt)
            mean_deltas   = mid_values(deltas)
            #print(deltas)
            W_delta       = weight * np.reciprocal(mean_deltas)

            
            
            first_climb = True
            first_accel = True


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
            thrust_type = 'HiTemp'

            
            if thrust_type == 'HiTemp':
                midpoint_Temps = temperature_profile(mid_values(alt), Tflex)
            else:
                midpoint_Temps = temperature_profile(mid_values(alt), Tapt)

            sins_gamma = np.zeros(np.shape(steps))
            segment_accel = np.zeros(np.shape(steps))
            ROCs = np.multiply(np.diff(alt), np.reciprocal(np.diff(times))) * 60.
            self.ROCs = ROCs #ft/min

            for i, step in enumerate(steps):
                if step == 'TakeOff':
                    

                    B8 = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            TO_Flap_ID[0]]['B'].iloc[0]
                    theta = temperature_ratio_(Alt=0, Tapt=Tapt)
                    
                    W_delta_TO = weight / pressure_ratio_(Alt=0, Papt=Papt)
                    
                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                            thrust_rating[thrust_type]['TO']].reset_index(drop=True)
                    
                    if thrust_type == 'HiTemp':
                        T_thrust = Tflex
                    else:
                        T_thrust = Tapt

                    #Should be CAS of first point
                    Fn_delta = corrected_net_thrust(thrust_coeff,mid_values(cas)[i+1],
                        temperature_profile(Alt=0, Tapt=T_thrust), Alt=0, Papt=Papt)
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
                    radar_first_CAS = mid_values(cas)[i]


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
                elif step == 'Accelerate' and first_accel:
                    

                    R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            accel_Flap_ID[0]]['R'].iloc[0]

                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['C']].reset_index(drop=True)

                    Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                    segment_accel[i] = distance_accelerate(tas_diff[i],
                        tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])

                    first_accel = False
                    second_accel = True
                elif step == 'Accelerate' and first_accel == False and \
                        second_accel==True:

                    R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            accel_Flap_ID[1]]['R'].iloc[0]

                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['C']].reset_index(drop=True)

                    Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                    segment_accel[i] = distance_accelerate(tas_diff[i],
                        tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
                 
                    second_accel = False
                elif step == 'Accelerate' and first_accel == False and \
                        second_accel==False:

                    R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            accel_Flap_ID[2]]['R'].iloc[0]

                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['C']].reset_index(drop=True)

                    Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                    segment_accel[i] = distance_accelerate(tas_diff[i],
                        tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
                    
        
            sins_gamma_est = sins_gamma_estimation(dist, alt, cas)
            
            
            #sins_gamma_o = (sins_gamma - min(sins_gamma[steps=='Climb']))/\
            #        (max(sins_gamma[steps=='Climb'])-min(sins_gamma[steps=='Climb']))
            #sins_gamma_est_o = (sins_gamma_est - min(sins_gamma_est[steps=='Climb']))/\
            #        (max(sins_gamma_est[steps=='Climb'])-min(sins_gamma_est[steps=='Climb']))
            #
            #seg_accel_o = (segment_accel[steps=='Accelerate']- \
            #        min(segment_accel[steps=='Accelerate'])) / \
            #        (max(segment_accel[steps=='Accelerate']) - \
            #        min(segment_accel[steps=='Accelerate']))
            #seg_accel_RADAR_o = (np.diff(dist)[steps=='Accelerate'] - \
            #        min(np.diff(dist)[steps=='Accelerate']))        / \
            #        (max(np.diff(dist)[steps=='Accelerate'])        - \
            #        min(np.diff(dist)[steps=='Accelerate']))

            sins_gamma_o = sins_gamma[steps=='Climb']
            sins_gamma_est_o =sins_gamma_est[steps=='Climb']            
            seg_accel_o = segment_accel[steps=='Accelerate']
            seg_accel_RADAR_o = np.diff(dist)[steps=='Accelerate']
           
            
            #cost_climb = (sins_gamma_o-sins_gamma_est_o)**2

            #cost_climb = cost_climb / max(cost_climb)
            #        
            #cost_climb = np.sum(cost_climb)
            #cost_accel = (seg_accel_RADAR_o - seg_accel_o)**2

            #cost_accel = cost_accel / max(cost_accel) 
            #cost_accel = np.sum(cost_accel)

            cost_climb = np.mean((sins_gamma[steps=='Climb'] - \
                   sins_gamma_est[steps=='Climb'])**2)
            
            cost_accel = np.mean((segment_accel[steps=='Accelerate'] - \
                   np.diff(dist)[steps=='Accelerate'])**2)
            
            cost_first_CAS = (est_first_climb_CAS - radar_first_CAS)**2

            cost_TO8 = (est_TO8 - np.diff(dist)[steps=='TakeOff'])**2
            
            #cost = cost_climb + cost_accel + 2*cost_first_CAS/200 + cost_TO8[0]/1e8
            cost = cost_first_CAS *  cost_TO8 * cost_accel **2
            
            #cost = cost_accel

            costs.append(cost) 

            #print('ANP climb angle: '  + str(sins_gamma[steps=='Climb']))
            #print('RADAR climb angle ' +  str(sins_gamma_est[steps=='Climb']))
            #print('ANP accel segment '+ str(segment_accel[steps=='Accelerate']))
            #print('RADAR accel segment ' + str(np.diff(dist)[steps=='Accelerate']))
            #print('ANP first CAS ' + str(est_first_climb_CAS))
            #print('RADAR first CAS ' + str(radar_first_CAS))
            #print('ANP TO8 ' + str(est_TO8))
            #print('RADAR TO8 ' + str(np.diff(dist)[0]))
            return cost 
        
                
        stage_length = 1
        weight = self.Default_weights[self.Default_weights['Stage Length']==\
        stage_length]['Weight (lb)'].iloc[0]
        Tflex = 60. #Celsius 

        x0 = [weight, Tflex]
        
        bounds =  [(weight - .2*weight, weight + .5*weight),
                   (Tflex -  .5*Tflex , Tflex  + .0*Tflex)]
        new_profile = optimize.differential_evolution(loss, bounds)
        new_profile2 = optimize.minimize(loss, x0, bounds=bounds)
        self.weight_GA = new_profile.x[0]
        self.Tflex_GA  = new_profile.x[1]
        self.weight_grad = new_profile2.x[0]
        self.Tflex_grad  = new_profile2.x[1]
                #sins_gamma_est = sins_gamma_estimation(dist, alt, cas)

                #print(steps)
                #print('ANP climb angle: '  + str(sins_gamma[steps=='Climb']))
                #print('RADAR climb angle ' +  str(sins_gamma_est[steps=='Climb']))
                #print('ANP accel segment '+ str(segment_accel[steps=='Accelerate']))
                #print('RADAR accel segment ' + str(np.diff(dist)[steps=='Accelerate']))
                #print('ANP first CAS ' + str(est_first_climb_CAS))
                #print('RADAR first CAS ' + str(radar_first_CAS))
                #print('ANP TO8 ' + str(est_TO8))
                #print('RADAR TO8 ' + str(np.diff(dist)[0]))

    #def segment_time(self, dist, dist_full, times):
    #    ptimes = np.interp(dist, dist_full, times)
    #    return ptimes

    def new_vert_profile_pinv(self, column_names, op_type = 'D'):
        
        print('Laying down new profile...')
        costs = []
        
        #weight = 127496.7607

        dist   = self.d_seg
        alt    = self.h_seg
        cas    = self.cas_seg
        times  = self.time_seg
        steps = self.steps
        #print('Laying down new profile...')
        
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

        #Tflex = (Tflex * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)

        
        Aero_coeff = self.Aerodynamic_coefficients
        dep_Aero_coeff = Aero_coeff[Aero_coeff["Op Type"] == op_type]
        Tapt = (Tapt * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)
        Papt = Papt / 33.864 #Conversion mbar to inHg (not ideal)
        

        sigmas       = air_density_ratio(alt, Tapt, Papt)
        tas          = np.multiply(cas, np.reciprocal(np.sqrt(sigmas)))
        mean_sigmas  = mid_values(sigmas)
        mean_sigmas_ = air_density_ratio(mid_values(alt), Tapt, Papt)

        tas_geom_mean = np.sqrt(mid_values(np.power(tas, 2))) #Impact eq-17
        tas_diff      = np.diff(np.power(tas, 2))
        deltas        = pressure_ratio_(alt, Papt)
        mean_deltas   = mid_values(deltas)
        #print(deltas)
        #W_delta       = weight * np.reciprocal(mean_deltas)

        
        
        first_climb = True
        first_accel = True


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
        thrust_type = 'HiTemp'

        
        #if thrust_type == 'HiTemp':
        #    midpoint_Temps = temperature_profile(mid_values(alt), Tflex)
        #else:
        #    midpoint_Temps = temperature_profile(mid_values(alt), Tapt)

        sins_gamma = np.zeros(np.shape(steps))
        segment_accel = np.zeros(np.shape(steps))
        ROCs = np.multiply(np.diff(alt), np.reciprocal(np.diff(times))) * 60.
        self.ROCs = ROCs #ft/min
        
        sins_gamma_RADAR = sins_gamma_estimation(dist, alt, cas)
        accel_seg_RADAR  = np.diff(dist)
        
        A = np.zeros((len(steps) - 1 + 1, 2)) #-1 no takeoff +1 first CAS
        b = np.zeros((len(steps) - 1 + 1, 1)) #-1 no takeoff +1 first CAS

        for i, step in enumerate(steps):
            if step == 'TakeOff':
                
                pass
                #B8 = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                #        TO_Flap_ID[0]]['B'].iloc[0]
                #theta = temperature_ratio_(Alt=0, Tapt=Tapt)
                #
                #W_delta_TO = weight / pressure_ratio_(Alt=0, Papt=Papt)
                #
                #thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                #        thrust_rating[thrust_type]['TO']].reset_index(drop=True)
                #
                #if thrust_type == 'HiTemp':
                #    T_thrust = Tflex
                #else:
                #    T_thrust = Tapt

                ##Should be CAS of first point
                #Fn_delta = corrected_net_thrust(thrust_coeff,mid_values(cas)[i+1],
                #    temperature_profile(Alt=0, Tapt=T_thrust), Alt=0, Papt=Papt)
                #est_TO8 = TO_ground_roll_distance(B8, theta, W_delta_TO, N, 
                #        Fn_delta)
            
            elif step == 'Climb' and first_climb:
                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        TO_Flap_ID[0]]['R'].iloc[0]

                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['TO']].reset_index(drop=True)

                #Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                #        midpoint_Temps[i], mid_values(alt)[i], Papt)
                ##print(Fn_delta)
                #sins_gamma[i] = sin_of_climb_angle(N, Fn_delta, W_delta[i], R,
                #        mid_values(cas)[i])

                alpha = sin_of_climb_angle_pinv(N, sins_gamma_RADAR[i], R,
                        mid_values(cas)[i])

                beta, H = corrected_net_thrust_pinv(thrust_coeff, 
                        mid_values(cas)[i], mid_values(alt)[i], Papt)

                A[i-1, :] = [alpha / mean_deltas[i], -H]
                b[i-1]      = beta

                #FIRST CAS

                C = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        TO_Flap_ID[0]]['C'].iloc[0]

                #est_first_climb_CAS = first_climb_CAS(C, weight)
                radar_first_CAS = mid_values(cas)[i]
                
                ksi = first_climb_CAS_pinv(C, radar_first_CAS)
                
                A[-1, :] = [1, 0]
                b[-1]      = ksi
                
                first_climb = False


            elif step == 'Climb' and first_climb==False:

                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        climb_Flap_ID[1]]['R'].iloc[0]
                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['C']].reset_index(drop=True)

               # Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
               #         midpoint_Temps[i], mid_values(alt)[i], Papt)

               # sins_gamma[i] = sin_of_climb_angle(N, Fn_delta, W_delta[i], R,
               #         mid_values(cas)[i])
                
                alpha = sin_of_climb_angle_pinv(N, sins_gamma_RADAR[i], R,
                        mid_values(cas)[i])

                beta, H = corrected_net_thrust_pinv(thrust_coeff, 
                        mid_values(cas)[i], mid_values(alt)[i], Papt)

                A[i-1, :] = [alpha / mean_deltas[i], -H]
                b[i-1]      = beta
            
            elif step == 'Accelerate' and first_accel:
                

                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        accel_Flap_ID[0]]['R'].iloc[0]

                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                    thrust_rating[thrust_type]['C']].reset_index(drop=True)

               # Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
               #     midpoint_Temps[i], mid_values(alt)[i], Papt)
               # segment_accel[i] = distance_accelerate(tas_diff[i],
               #     tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])

                phi = distance_accel_pinv(tas_diff[i], tas_geom_mean[i], N, 
                        accel_seg_RADAR[i], R, ROCs[i])

                beta, H = corrected_net_thrust_pinv(thrust_coeff, 
                        mid_values(cas)[i], mid_values(alt)[i], Papt)

                A[i-1, :] = [phi / mean_deltas[i], -H]
                b[i-1]      = beta
                
                
                first_accel = False
                second_accel = True
            elif step == 'Accelerate' and first_accel == False and \
                    second_accel==True:

                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        accel_Flap_ID[1]]['R'].iloc[0]

                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                    thrust_rating[thrust_type]['C']].reset_index(drop=True)

                #Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                #    midpoint_Temps[i], mid_values(alt)[i], Papt)
                #segment_accel[i] = distance_accelerate(tas_diff[i],
                #    tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
                
                phi = distance_accel_pinv(tas_diff[i], tas_geom_mean[i], N, 
                        accel_seg_RADAR[i], R, ROCs[i])

                beta, H = corrected_net_thrust_pinv(thrust_coeff, 
                        mid_values(cas)[i], mid_values(alt)[i], Papt)

                A[i-1, :] = [phi / mean_deltas[i], -H]
                b[i-1]      = beta
             
                second_accel = False
            elif step == 'Accelerate' and first_accel == False and \
                    second_accel==False:

                R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                        accel_Flap_ID[2]]['R'].iloc[0]

                thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                    thrust_rating[thrust_type]['C']].reset_index(drop=True)

                #Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                #    midpoint_Temps[i], mid_values(alt)[i], Papt)
                #segment_accel[i] = distance_accelerate(tas_diff[i],
                #    tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
                    
                phi = distance_accel_pinv(tas_diff[i], tas_geom_mean[i], N, 
                        accel_seg_RADAR[i], R, ROCs[i])

                beta, H = corrected_net_thrust_pinv(thrust_coeff, 
                        mid_values(cas)[i], mid_values(alt)[i], Papt)

                A[i-1, :] = [phi / mean_deltas[i], -H]
                b[i-1]      = beta
        
            
        x = np.dot(np.linalg.pinv(A), b)       
        self.weight_pinv = x[0]
        self.Tflex_pinv  = x[1]
    def plot_ANP_profile(self, column_names, op_type = 'D'):
        def model(x):
            weight = x[0]
            Tflex  = x[1]


            dist   = self.d_seg
            alt    = self.h_seg
            cas    = self.cas_seg
            times  = self.time_seg
            steps = self.steps
            
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

            Tflex = (Tflex * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)

            
            Aero_coeff = self.Aerodynamic_coefficients
            dep_Aero_coeff = Aero_coeff[Aero_coeff["Op Type"] == op_type]
            Tapt = (Tapt * 9/5) + 32.0 #Conversion to Farenheit from C (not ideal)
            Papt = Papt / 33.864 #Conversion mbar to inHg (not ideal)
            

            sigmas       = air_density_ratio(alt, Tapt, Papt)
            tas          = np.multiply(cas, np.reciprocal(np.sqrt(sigmas)))
            mean_sigmas  = mid_values(sigmas)
            mean_sigmas_ = air_density_ratio(mid_values(alt), Tapt, Papt)

            tas_geom_mean = np.sqrt(mid_values(np.power(tas, 2))) #Impact eq-17
            tas_diff      = np.diff(np.power(tas, 2))
            deltas        = pressure_ratio_(alt, Papt)
            mean_deltas   = mid_values(deltas)
            #print(deltas)
            W_delta       = weight * np.reciprocal(mean_deltas)

            
            
            first_climb = True
            first_accel = True


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
            thrust_type = 'HiTemp'

            
            if thrust_type == 'HiTemp':
                midpoint_Temps = temperature_profile(mid_values(alt), Tflex)
            else:
                midpoint_Temps = temperature_profile(mid_values(alt), Tapt)

            sins_gamma = np.zeros(np.shape(steps))
            segment_accel = np.zeros(np.shape(steps))
            ROCs = np.multiply(np.diff(alt), np.reciprocal(np.diff(times))) * 60.
            self.ROCs = ROCs #ft/min

            for i, step in enumerate(steps):
                if step == 'TakeOff':
                    

                    B8 = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            TO_Flap_ID[0]]['B'].iloc[0]
                    theta = temperature_ratio_(Alt=0, Tapt=Tapt)
                    
                    W_delta_TO = weight / pressure_ratio_(Alt=0, Papt=Papt)
                    
                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                            thrust_rating[thrust_type]['TO']].reset_index(drop=True)
                    
                    if thrust_type == 'HiTemp':
                        T_thrust = Tflex
                    else:
                        T_thrust = Tapt

                    #Should be CAS of first point
                    Fn_delta = corrected_net_thrust(thrust_coeff,mid_values(cas)[i+1],
                        temperature_profile(Alt=0, Tapt=T_thrust), Alt=0, Papt=Papt)
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
                elif step == 'Accelerate' and first_accel:
                    

                    R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            accel_Flap_ID[0]]['R'].iloc[0]

                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['TO']].reset_index(drop=True)

                    Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                    segment_accel[i] = distance_accelerate(tas_diff[i],
                        tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])

                    first_accel = False
                    second_accel = True
                elif step == 'Accelerate' and first_accel == False and \
                        second_accel==True:

                    R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            accel_Flap_ID[1]]['R'].iloc[0]

                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['TO']].reset_index(drop=True)

                    Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                    segment_accel[i] = distance_accelerate(tas_diff[i],
                        tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
                 
                    second_accel = False
                elif step == 'Accelerate' and first_accel == False and \
                        second_accel==False:

                    R = dep_Aero_coeff[dep_Aero_coeff['Flap_ID']==\
                            accel_Flap_ID[2]]['R'].iloc[0]

                    thrust_coeff = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']==\
                        thrust_rating[thrust_type]['TO']].reset_index(drop=True)

                    Fn_delta = corrected_net_thrust(thrust_coeff, mid_values(cas)[i],
                        midpoint_Temps[i], mid_values(alt)[i], Papt)
                    segment_accel[i] = distance_accelerate(tas_diff[i],
                        tas_geom_mean[i], N, Fn_delta, W_delta[i], R, ROCs[i])
                    
        
            sins_gamma_est = sins_gamma_estimation(dist, alt, cas)
            
            
            sins_gamma_est_o =sins_gamma_est[steps=='Climb']            
            seg_accel_o = segment_accel[steps=='Accelerate']

            return sins_gamma_est, segment_accel, est_first_climb_CAS, est_TO8,\
                    cas

        sins, accels, first_cas, to8, cas = model([self.weight_GA, self.Tflex_GA])
        sins_g, accels_g, first_cas_g, to8_g, cas_g =\
                model([self.weight_grad, self.Tflex_grad])
        
        sins_p, accels_p, first_cas_p, to8_p, cas_p =\
                model([self.weight_pinv, self.Tflex_pinv])
        
        
        weight = ((self.cas_seg[1:] + self.cas_seg[:-1])[1]/2 / 0.394884)**2
        Tflex = 50
        sins_n, accels_n, first_cas_n, to8_n, cas_n =\
                model([weight, Tflex])
        
        


        for i, step in enumerate(self.steps):
            if step == 'TakeOff':
                d_anp = np.r_[0, to8]
                h_anp = np.r_[0, 0]
                cas_anp = np.r_[0, first_cas]
            elif step == 'Climb':
                h_anp = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp = np.r_[d_anp, d_anp[-1] + (np.diff(self.h_seg)[i] \
                        /np.tan(np.arcsin(sins[i])))]
                cas_anp = np.r_[cas_anp, cas_anp[-1]]
            elif step == 'Accelerate':
                h_anp   = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp   = np.r_[d_anp, d_anp[-1] + accels[i]]
                cas_anp = np.r_[cas_anp, cas_anp[-1] + np.diff(cas)[i]]
        
        
        self.d_anp = np.copy(d_anp)
        self.h_anp = np.copy(h_anp)
        self.cas_anp = np.copy(cas_anp)
        
        d_anp = []
        h_anp = []
        cas_anp = []
        for i, step in enumerate(self.steps):
            if step == 'TakeOff':
                d_anp = np.r_[0, to8_p]
                h_anp = np.r_[0, 0]
                cas_anp = np.r_[0, first_cas_p]
            elif step == 'Climb':
                h_anp = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp = np.r_[d_anp, d_anp[-1] + (np.diff(self.h_seg)[i] \
                        /np.tan(np.arcsin(sins_p[i])))]
                cas_anp = np.r_[cas_anp, cas_anp[-1]]
            elif step == 'Accelerate':
                h_anp   = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp   = np.r_[d_anp, d_anp[-1] + accels_p[i]]
                cas_anp = np.r_[cas_anp, cas_anp[-1] + np.diff(cas_p)[i]]
        
        self.d_anp_p = np.copy(d_anp)
        self.h_anp_p = np.copy(h_anp)
        self.cas_anp_p = np.copy(cas_anp)
        
        d_anp = []
        h_anp = []
        cas_anp = []
        
        for i, step in enumerate(self.steps):
            if step == 'TakeOff':
                d_anp = np.r_[0, to8_n]
                h_anp = np.r_[0, 0]
                cas_anp = np.r_[0, first_cas_n]
            elif step == 'Climb':
                h_anp = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp = np.r_[d_anp, d_anp[-1] + (np.diff(self.h_seg)[i] \
                        /np.tan(np.arcsin(sins_n[i])))]
                cas_anp = np.r_[cas_anp, cas_anp[-1]]
            elif step == 'Accelerate':
                h_anp   = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp   = np.r_[d_anp, d_anp[-1] + accels_n[i]]
                cas_anp = np.r_[cas_anp, cas_anp[-1] + np.diff(cas_p)[i]]
        
        self.d_anp_n = np.copy(d_anp)
        self.h_anp_n = np.copy(h_anp)
        self.cas_anp_n = np.copy(cas_anp)
        
        d_anp = []
        h_anp = []
        cas_anp = []
        for i, step in enumerate(self.steps):
            if step == 'TakeOff':
                d_anp = np.r_[0, to8_g]
                h_anp = np.r_[0, 0]
                cas_anp = np.r_[0, first_cas_g]
            elif step == 'Climb':
                h_anp = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp = np.r_[d_anp, d_anp[-1] + (np.diff(self.h_seg)[i] \
                        /np.tan(np.arcsin(sins_g[i])))]
                cas_anp = np.r_[cas_anp, cas_anp[-1]]
            elif step == 'Accelerate':
                h_anp   = np.r_[h_anp, h_anp[-1] + np.diff(self.h_seg)[i]]
                d_anp   = np.r_[d_anp, d_anp[-1] + accels_g[i]]
                cas_anp = np.r_[cas_anp, cas_anp[-1] + np.diff(cas_p)[i]]
        
        self.d_anp_g = np.copy(d_anp)
        self.h_anp_g = np.copy(h_anp)
        self.cas_anp_g = np.copy(cas_anp)
        
        def plot_anp():
            px = self.d_seg
            py = self.h_seg
            pcas = self.cas_seg
            X = self.d
            Y = self.h
            CAS = self.cas
            
            plt.subplot(211)
            plt.plot(X, Y, '.', label='Radar')
            #plt.plot(px, py, '-or', label='Model')
            plt.plot(self.d_anp, self.h_anp, '-og', label='GA')
            plt.plot(self.d_anp_p, self.h_anp_p, '-oy', label='pinv')
            plt.plot(self.d_anp_g, self.h_anp_g, '-om', label='grad')
            plt.xlabel('Dist [ft]')
            plt.ylabel('Alt [ft]')
            plt.legend()
            
            plt.subplot(212)
            plt.plot(X, CAS, '.', label='Radar')
            #plt.plot(px, pcas, '-or', label='Model')
            plt.plot(self.d_anp, self.cas_anp, '-og', label='GA')
            plt.plot(self.d_anp_p, self.cas_anp_p, '-oy', label='pinv')
            plt.plot(self.d_anp_g, self.cas_anp_g, '-om', label='grad')
            plt.xlabel('Dist [ft]')
            plt.ylabel('CAS [kts]')
            #plt.legend()

            
            plt.show()

        plot_anp()
