import pandas as pd 
import numpy as np
from scipy import optimize, interpolate
import matplotlib.pyplot as plt
import pdb
from sklearn.isotonic import IsotonicRegression
def get_data():
    flight1 = pd.read_csv(r'data/radar_tracks/ORY_to_BIA/flight_augmented.csv')
    #flight1 = flights[flights['Flight id'] == 1]
    xdata = flight1['Distance (ft)'].to_numpy()
    ydata = flight1['altitude (ft)'].to_numpy()
    cas = flight1['CAS (kts)'].to_numpy()
    xdata = xdata[22:70]
    ydata = np.float64(ydata[22:70])
    cas = cas[22:70]
    cas_o = cas[ydata>0]
    xdata_o = xdata[ydata>0]
    ydata_o = ydata[ydata>0]
    return xdata, ydata, cas, cas_o, xdata_o, ydata_o

def segments_fit(X, Y, cas, maxcount):
    xmin = X.min()
    xmax = X.max()

    n = len(X)

    AIC_ = float('inf')
    BIC_ = float('inf')
    r_ = None
    
    for count in [4]:#range(1, maxcount + 1):
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

       # AIC = n * np.log10(err(r.x)) + 4 * count
       # BIC = n * np.log10(err(r.x)) + 2 * count * np.log(n)
       # 
       # if ((BIC < BIC_) & (AIC < AIC_)):
       #     r_ = r
       #     AIC_ = AIC
       #     BIC_ = BIC
       # else: 
       #     count = count - 1
       #     break
    
    return func(r.x), r

X, Y, CAS, cas_o, x_o, y_o = get_data()
vals, res = segments_fit(x_o, y_o, cas_o, 8)
px, py, pcas = vals
plt.subplot(211)
plt.plot(X, Y, '.')
plt.plot(px, py, '-or')
plt.subplot(212)
plt.plot(X, CAS, '.')
plt.plot(px, pcas, '-or')
plt.show()
