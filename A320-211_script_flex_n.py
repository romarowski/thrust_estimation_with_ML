import pandas as pd 
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import pdb
def get_data():
    flights = pd.read_csv(r'data/radar_tracks/A320_211.csv')
    flight1 = flights[flights['Flight id'] == 1]
    xdata = flight1['Distance (ft)'].to_numpy()
    ydata = flight1['Altitude (ft)'].to_numpy()
    xdata = xdata[2:]
    ydata = ydata[2:]
    return xdata, ydata

def segments_fit(X, Y, maxcount):
    #pdb.set_trace()
    xmin = X.min()
    xmax = X.max()

    n = len(X)

    AIC_ = float('inf')
    BIC_ = float('inf')
    r_ = None
    
    for count in range(1, maxcount + 1):
        seg = np.full(count - 1, (xmax - xmin) / count)

        px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        mask = [[np.abs(X - x) < (xmax - xmin) * 0.1] for x in px_init]
        py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.1].mean() for x in px_init])

        def func(p):
            pdb.set_trace()
            seg = p[:count - 1]
            py = p[count - 1:]
            px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
            return px, py

        def err(p):
            px, py = func(p)
            Y2 = np.interp(X, px, py)
            return np.mean((Y - Y2)**2)

        r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')

        AIC = n * np.log10(err(r.x)) + 4 * count
        BIC = n * np.log10(err(r.x)) + 2 * count * np.log(n)
        
        if ((BIC < BIC_) & (AIC < AIC_)):
            r_ = r
            AIC_ = AIC
            BIC_ = BIC
        else: 
            count = count - 1
            break
    
    return func(r_.x)

X, Y = get_data()
px, py = segments_fit(X, Y, 8)

plt.plot(X, Y, '.')
plt.plot(px, py, '-or')
plt.show()
