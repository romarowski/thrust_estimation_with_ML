import time
import pandas as pd 
import numpy as np
from scipy import optimize, interpolate
import matplotlib.pyplot as plt
import pdb
from sklearn.isotonic import IsotonicRegression
from modules import NoiseCraft
def get_data():
    flight1 = pd.read_csv(r'data/radar_tracks/ORY_to_BIA/flight_augmented.csv')
    #flight1 = flights[flights['Flight id'] == 1]
    xdata = flight1['Distance (ft)'].to_numpy()
    ydata = flight1['altitude (ft)'].to_numpy()
    cas = flight1['CAS (kts)'].to_numpy()
    time = flight1['snapshot'].to_numpy()
    time = time[43:71]
    xdata = xdata[43:71]
    ydata = np.float64(ydata[43:71])
    cas = cas[43:71]
    cas_o = cas[ydata>0]
    xdata_o = xdata[ydata>0]
    ydata_o = ydata[ydata>0]
    return xdata, ydata, cas, cas_o, xdata_o, ydata_o, time

def segments_fit(X, Y, cas, maxcount):
    xmin = X.min()
    xmax = X.max()

    n = len(X)

    AIC_ = float('inf')
    BIC_ = float('inf')
    r_ = None
    
    for count in [4]: #in range(1, maxcount + 1): #[1, 2,  4, 6]:
        seg = np.full(count - 1, (xmax - xmin) / count)

        px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        mask = [[np.abs(X - x) < (xmax - xmin) * 0.1] for x in px_init]
         

        py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.1].mean() for x in px_init])
        
        pcas_init = np.array([cas[np.abs(X - x) < (xmax - xmin) * 0.1].mean() for x in px_init])

        def func(p):
            seg = p[:count - 1]
            py = p[count - 1:2*count]
            #pcas = p[2*count:-2]
            pcas = p[2*count:]
            px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
            return px, py, pcas

        def err(p):
            px, py, pcas = func(p)
            #weights = p[-2:]
            #weights = weights/np.linalg.norm(weights)
            Y2 = np.interp(X, px, py)
            CAS2 = np.interp(X, px, pcas)
            if not np.all(np.diff(pcas)>=0):
                penalty = float('inf')
            else:
                penalty = 0
            #ir = IsotonicRegression(out_of_bounds='clip')
            #CAS2_ = ir.fit_transform(px, pcas)
            #CAS2 = ir.predict(X)
            #cas_interpolant = interpolate.PchipInterpolator(px, pcas)
            #CAS2 = cas_interpolant.__call__(X)
            #return weights[0]*np.mean((Y - Y2)**2) +\
            #       weights[1]*np.mean((cas - CAS2)**2)
            return np.mean((Y - Y2)**2) +\
                   1000*np.mean((cas - CAS2)**2) + penalty

        
        #x0=np.r_[seg, py_init, pcas_init, [.5,.5]]
        x0=np.r_[seg, py_init, pcas_init]
        r = optimize.minimize(err, x0=x0, method='Nelder-Mead')

        AIC = n * np.log10(err(r.x)) + 4 * count
        BIC = n * np.log10(err(r.x)) + 2 * count * np.log(n)
        vals = func(r.x)
        if AIC < AIC_: #(BIC < BIC_) and (AIC < AIC_):
            r_ = r
            AIC_ = AIC
            BIC_ = BIC
        #else: 
        #    count = count - 1
            #break
        #X_plot, Y_plot, CAS_plot, cas_o, x_o, y_o = get_data()
        #px, py, pcas = vals
        #
        #plt.subplot(211)
        #plt.plot(X_plot, Y_plot, '.')
        #plt.plot(px, py, '-or')
        #plt.ylabel('Altitude [ft]')
        #plt.subplot(212)
        #plt.xlabel('Distance [ft]')
        #plt.ylabel('CAS [kts]')
        #plt.plot(X_plot, CAS_plot, '.')
        #plt.plot(px, pcas, '-or')
        #plt.show()
    return func(r_.x), r_

def procedural_steps_takeoff(segmented_input, acraft, treshold_cas=10):
    dist = segmented_input[0]
    alt  = segmented_input[1]
    cas  = segmented_input[2]
    cas_derivative = np.diff(cas) <= treshold_cas
    steps = ['Climb' if der else 'Accelerate' for der in cas_derivative]
    climbs_deltaH = np.diff(alt)[cas_derivative]
    climbs_deltaS  = np.diff(dist)[cas_derivative]
    sins_gamma = climbs_deltaH / np.sqrt(climbs_deltaH**2 + climbs_deltaS**2)
    return steps, sins_gamma, cas[0]

def climb_at_ct_speed(x, params, acraft):
    
    Fn = x[0]
    W  = x[1]
    
    # CLIMBS WITH EQUATION B-12
    N = 2
    R = 0.071403
    #R = np.array([0.066827, 0.071403, 0.056281])
    K = 0.95
    bank_angle = 0
    sin_gamma = K * (N * Fn / W - R / np.cos(bank_angle)) 
    sin_gamma_radar = params['sins_gamma']
    cost = np.mean((sin_gamma - sin_gamma_radar) ** 2)
    
    #INITIAL CLIMB WITH EQ B-15
    C = TO_Aerodyn_coeff(acraft)
    weight = (params['first_CAS']/C)**2
    cost_first_climb = (W - weight)**2
    return cost + cost_first_climb 

def TO_Aerodyn_coeff(acraft, TO_Flap_ID = '1+F'):
    Aerodyn_coeff = acraft.Aerodynamic_coefficients
    C = Aerodyn_coeff[Aerodyn_coeff['Flap_ID']==TO_Flap_ID]['C'].iloc[0]
    return C

def moving_mean(vector):
    return (vector[:-1] + vector[1:]) / 2


X, Y, CAS, cas_o, x_o, y_o, times = get_data()
times -= times[0]
t = time.time()
vals, res = segments_fit(x_o, y_o, cas_o, 6)
elapsed = time.time() - t
print(elapsed)
noisecraft = NoiseCraft.NoiseCraft('A320-211')
steps, sins_gamma, first_CAS = procedural_steps_takeoff(vals, noisecraft)
params = {'sins_gamma' : sins_gamma,
          'first_CAS'  : first_CAS}
weight0 = 139200 #lb
Fn0 = 23255 #lb 
x0 = [Fn0, weight0]
#Jet_eng_coeff = noisecraft.Jet_engine_coefficients
#Max_climb_coeffs = Jet_eng_coeff[Jet_eng_coeff['Thrust Rating']=='Max Takeoff']
#E = Max_climb_coeffs['E'].iloc[0]
#F =
residual_first_climb = optimize.minimize(climb_at_ct_speed, x0, 
        args=(params, noisecraft))
px, py, pcas = vals
ptimes = np.interp(px, X, times)
plt.subplot(221)
plt.plot(X, Y, '.')
plt.plot(px, py, '-or')
plt.subplot(223)
plt.plot(X, CAS, '.')
plt.plot(px, pcas, '-or')
plt.subplot(222)
plt.plot(times, Y, '.')
plt.plot(ptimes, py, '-or')
plt.subplot(224)
plt.plot(times, CAS, '.')
plt.plot(ptimes, pcas, '-or')
plt.show()
