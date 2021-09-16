import pandas as pd
import plotly.express as px
import numpy as np
from scipy.optimize import curve_fit, minimize, shgo
import matplotlib.pyplot as plt
from matplotlib import animation
import pdb
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML
def func_(d, a1,b1,mu1,a2,b2,mu2,a3,b3,mu3,a4,b4,mu4):
    return (a1*(d-mu1)+b1)*(np.heaviside(d-mu1,0)-np.heaviside(d-mu2,0))+\
           (a2*(d-mu2)+b2)*(np.heaviside(d-mu2,0)-np.heaviside(d-mu3,0))+\
           (a3*(d-mu3)+b3)*(np.heaviside(d-mu3,0)-np.heaviside(d-mu4,0))+\
           (a4*(d-mu4)+b4)*(np.heaviside(d-mu4,0))

def func(d, a1,b1,mu1,a2,b2,mu2,a3,b3,mu3,a4,b4,mu4):
    y = np.piecewise(d,
            [d<mu1,d>=mu1 and d<mu2,d>=mu2 and d<mu3,d>=mu3],
            [lambda d:a1*(d-mu1)+b1,
             lambda d:a2*(d-mu2)+b2,
             lambda d:a3*(d-mu3)+b3,
             lambda d:a4*(d-mu4)+b4
             ]
            )
    return y
def residuals(x, d, altitude, cas):
    a1  = x[0]
    mu1 = x[1] 
    b1  = x[2]
    a2  = x[3]
    mu2 = x[4]
    b2  = x[5]
    a3  = x[6]
    mu3 = x[7]
    b3  = x[8]
    a4  = x[9]
    mu4 = x[10]
    b4  = x[11]
    c1 = x[12]
    d1 = x[13]
    c2 = x[14]
    d2 = x[15]
    c3 = x[16]
    d3 = x[17]
    c4 = x[18]
    d4 = x[19]
    alpha = x[20]
    beta = x[21]
    alpha = 1
    beta =  1000
    #mu1=0
    alt_hat =\
            (a1*(d-mu1)+b1)*(np.heaviside(d-mu1,0)-np.heaviside(d-mu2,0))+\
            (a2*(d-mu2)+b2)*(np.heaviside(d-mu2,0)-np.heaviside(d-mu3,0))+\
            (a3*(d-mu3)+b3)*(np.heaviside(d-mu3,0)-np.heaviside(d-mu4,0))+\
            (a4*(d-mu4)+b4)*(np.heaviside(d-mu4,0))
    cas_hat =\
    (c1*(d-mu1)+d1)*(np.heaviside(d-mu1,0)-np.heaviside(d-mu2,0))+\
    (c2*(d-mu2)+d2)*(np.heaviside(d-mu2,0)-np.heaviside(d-mu3,0))+\
    (c3*(d-mu3)+d3)*(np.heaviside(d-mu3,0)-np.heaviside(d-mu4,0))+\
    (c4*(d-mu4)+d4)*(np.heaviside(d-mu4,0))
 
    
    return alpha*np.mean((alt_hat - altitude)**2) +\
           beta *np.mean((cas - cas_hat)**2) 


def model(x, d):
    a1  = x[0]
    mu1 = x[1] 
    b1  = x[2]
    a2  = x[3]
    mu2 = x[4]
    b2  = x[5]
    a3  = x[6]
    mu3 = x[7]
    b3  = x[8]
    a4  = x[9]
    mu4 = x[10]
    b4  = x[11]
    c1 = x[12]
    d1 = x[13]
    c2 = x[14]
    d2 = x[15]
    c3 = x[16]
    d3 = x[17]
    c4 = x[18]
    d4 = x[19]
    alpha = x[20]
    beta = x[21]
    alpha = 1
    beta=1
    #mu1=0
    alt_hat = (a1*(d-mu1)+b1)*(np.heaviside(d-mu1,0)-np.heaviside(d-mu2,0))+\
            (a2*(d-mu2)+b2)*(np.heaviside(d-mu2,0)-np.heaviside(d-mu3,0))+\
            (a3*(d-mu3)+b3)*(np.heaviside(d-mu3,0)-np.heaviside(d-mu4,0))+\
            (a4*(d-mu4)+b4)*(np.heaviside(d-mu4,0))
    cas_hat =\
    (c1*(d-mu1)+d1)*(np.heaviside(d-mu1,0)-np.heaviside(d-mu2,0))+\
            (c2*(d-mu2)+d2)*(np.heaviside(d-mu2,0)-np.heaviside(d-mu3,0))+\
            (c3*(d-mu3)+d3)*(np.heaviside(d-mu3,0)-np.heaviside(d-mu4,0))+\
            (c4*(d-mu4)+d4)*(np.heaviside(d-mu4,0))
 
    
    return alt_hat, cas_hat 

#def model(x, d):
#    a1  = x[0]
#    mu1 = x[1] 
#    b1  = x[2]
#    a2  = x[3]
#    mu2 = x[4]
#    b2  = x[5]
#    a3  = x[6]
#    mu3 = x[7]
#    b3  = x[8]
#    a4  = x[9]
#    mu4 = x[10]
#    b4  = x[11]
#    mu1=0
#    return (a1*(d-mu1)+b1)*(np.heaviside(d-mu1,0)-np.heaviside(d-mu2,0))+\
#           (a2*(d-mu2)+b2)*(np.heaviside(d-mu2,0)-np.heaviside(d-mu3,0))+\
#           (a3*(d-mu3)+b3)*(np.heaviside(d-mu3,0)-np.heaviside(d-mu4,0))+\
#           (a4*(d-mu4)+b4)*(np.heaviside(d-mu4,0))

def initialize(a1, a2, a3, a4, mu1, mu2, mu3, mu4, b1, 
               c1, c2, c3, c4, d1, alpha=0.5, beta=.5):
     b2 = a1*(mu2-mu1)+b1
     b3 = a2*(mu3-mu2)+b2
     b4 = a3*(mu4-mu3)+b3

     d2 = c1*(mu2-mu1)+d1
     d3 = c2*(mu3-mu2)+d2
     d4 = c3*(mu4-mu3)+d3
      
     x0 = [a1, mu1, b1, a2, mu2, b2, a3, mu3, b3, a4, mu4, b4,
           c1, d1, c2, d2, c3, d3, c4, d4, alpha, beta]
     return x0
 
def store(history):
    def store_it(xk):
        fname = r'./iterations_w_cas.txt'
        header = 'a1 mu1 b1 a2 mu2 b2 a3 mu3 b3 a4 mu4 b4'
        fmt = '%.2e'
        history.append(np.copy(xk))
        #np.savetxt(fname=fname, X=xk, fmt=fmt)
    return store_it
def get_data():
    flight1 = pd.read_csv(r'data/radar_tracks/ORY_to_BIA/flight_augmented.csv')
    xdata = flight1['Distance (ft)'].to_numpy()
    ydata = flight1['altitude (ft)'].to_numpy()
    cas = flight1['CAS (kts)'].to_numpy()
    xdata = xdata[22:70]
    ydata = np.float64(ydata[22:70])
    cas = cas[22:70]
    cas_o = cas[ydata>0]
    xdata_o = xdata[ydata>0]
    ydata_o = ydata[ydata>0]
        
    # TODO [NORMALIZE DATA]
    #pdb.set_trace()
    xmax   = np.max(xdata_o)
    ymax   = np.max(ydata_o)
    casmax = np.max(cas_o) 
    xdata_o /=np.max(xdata)
    ydata_o /= np.max(ydata)
    cas_o /= np.max(cas)
    return xdata, ydata, cas, xmax, ymax, casmax, cas_o, xdata_o, ydata_o



def main():
    mus = np.r_[0, np.random.rand(3)]
    mus.sort()
    x0 = initialize(a1=0, a2=np.random.rand(), a3=np.random.rand(), 
         a4=np.random.rand(),
         mu1=mus[0], mu2=mus[1], mu3=mus[2], mu4=mus[3],
         b1=0, 
         c1=np.random.rand(), c2=np.random.rand(), c3=np.random.rand(), 
         c4=np.random.rand(), d1=np.random.rand())
    eps = 0#1e-3
    cons = ({'type': 'ineq',
             'fun' : lambda x: np.array([x[1]+eps])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[4]+eps])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[7]+eps])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[10]+eps])},
          #  {'type': 'ineq',
          #   'fun' : lambda x: np.array([x[2]])},
          #  {'type': 'ineq',
          #   'fun' : lambda x: np.array([x[5]])},
          #  {'type': 'ineq',
          #   'fun' : lambda x: np.array([x[8]])},
          #  {'type': 'ineq',
          #   'fun' : lambda x: np.array([x[11]])},
            {'type': 'ineq',
             'fun': lambda x: np.array(x[:])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[4]+eps-x[1]])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[7]+eps-x[4]])},
            {'type': 'ineq',
            'fun' : lambda x: np.array([x[10]+eps-x[7]])},            
            {'type': 'eq',
             'fun' : lambda x: np.array([x[0]*(x[4]-x[1])+x[2]-x[5]])},
            {'type': 'eq',
             'fun' : lambda x: np.array([x[3]*(x[7]-x[4])+x[5]-x[8]])},
            {'type': 'eq',
             'fun' : lambda x: np.array([x[6]*(x[10]-x[7])+x[8]-x[11]])},
            {'type': 'eq',
             'fun' : lambda x: np.array([x[12]*(x[4]-x[1])+x[13]-x[15]])},
            {'type': 'eq',
             'fun' : lambda x: np.array([x[14]*(x[7]-x[4])+x[15]-x[17]])},
            {'type': 'eq',
             'fun' : lambda x: np.array([x[16]*(x[10]-x[7])+x[17]-x[19]])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[20]+eps])},
            {'type': 'ineq',
             'fun' : lambda x: np.array([x[21]+eps])},
             )
    #x0 = [.0, .0, .0, .1, 8000.0, .0, .2, 12e3, 500, .4, 15e3, 500]
    
    #Solution w 3 lines
    #x0 = np.array([ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
    #    4.20209911e-01, 3.89029345e+03,  0.00000000e+00,  
    #    4.61616055e-02,  2.95866930e+03, -3.91477700e+02, 
    #    2.20025059e-01,  1.71364203e+04,  2.62990047e+02])
    
    #Solution w 4 lines
   # x0 = np.array([ 2.71289513e-17,  1.36492607e+01, -7.78446012e-14,  1.85096795e-01,
   #     8.78525868e+03,  1.60119892e-13,  1.08964453e-01,  2.04519214e+04,
   #     2.15946188e+03,  5.35079596e-02,  2.77729771e+04,  2.95719671e+03])
    #x0 += np.random.randn(12,)
    xdata, ydata, cas, xmax, ymax, casmax, cas_o, x_o, y_o = get_data()
    history=[np.array(x0)]
    res = minimize(residuals, x0=x0, 
             args=(x_o, y_o, cas_o), constraints = cons, 
             callback=store(history))
    #bounds = [(0, .1), (0, .1),   (0, .1),
    #          (0.1, 2), (800, 1200), (0, max(ydata)),
    #          (.1, 2), (1200, max(xdata)), (0, max(ydata)),
    #          (.1, 3), (1200, max(xdata)), (0, max(ydata))] 
    
   #res = shgo(residuals, bounds = bounds, args=(xdata,ydata),
   #         constraints = cons) 
    #popt, pcov = curve_fit(func, xdata, ydata)

    alt_hat, cas_hat = model(res.x, x_o)
    x_o *= xmax
    y_o *= ymax
    cas_o   *= casmax
    alt_hat *= ymax
    cas_hat *= casmax
    fig_sol, ax1 = plt.subplots(2,1)
    ax1[0].plot(xdata, ydata, '-bo', label='Radar')
    ax1[0].plot(x_o, alt_hat, 'r', lw=2, label='Model')
    ax1[0].vlines(res.x[[1, 4, 7, 10]]*xmax, min(ydata), max(ydata), 
            linestyles='dashed', colors=['k', 'g', 'c', 'm'])
    ax1[0].set_xlabel('Distance')
    ax1[0].set_ylabel('Altitude')
    ax1[0].legend(loc='upper left')
    ax1[1].plot(xdata, cas, '-bo', label='Radar')
    ax1[1].plot(x_o, cas_hat, 'r', lw=2, label='Model')
    ax1[1].vlines(res.x[[1, 4, 7, 10]]*xmax, np.min(cas), np.max(cas), 
            linestyles='dashed', colors=['k', 'g', 'c', 'm'])
    ax1[1].set_xlabel('Distance')
    ax1[1].set_ylabel('CAS (kts)')
    ax1[1].legend(loc='upper left')
    plt.show()
   # count = 11
   # fig_sol.savefig('solution_{count}.png'.format(count=count), dpi=600)

   # fig, ax = plt.subplots()
   # ax.plot(xdata, ydata, 'bo-', label='Radar')
   # line, = ax.plot([], [], 'r', label='Model', lw=2)
   # ax.set_xlabel('Distance')
   # ax.set_ylabel('Altitude')
   ##history = np.loadtxt('iterations.txt') 
   # def init():
   #     line.set_data([], [])
   #     return line,
   # def animate(i, history, xdata):
   #     line.set_data(xdata, model(history[i], xdata))
   #     return line,
   # 
   # ax.set_xlim((min(ydata)-10, max(xdata)+10))
   # ax.set_ylim((min(ydata)-10, max(ydata)+10))

   # ax.legend(loc='upper left')
   # anim = animation.FuncAnimation(fig, animate, init_func=init, 
   #                                fargs = (history, xdata),
   #                                frames=res.nit, interval=1e3, blit=True,
   #                                repeat=False)
   # anim.save('opt_{count}.mp4'.format(count=count)
   #         , fps=1, extra_args=['-vcodec', 'libx264'])
    
    return res, xdata, ydata, x0, history
if __name__ == '__main__':
   res, xdata, ydata, x0, path = main()
   print(res.fun)
