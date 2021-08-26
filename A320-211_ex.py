from modules import corrected_net_thrust, CAS_TakeOff, TO_ground_roll_dist, \
        opt_TO_ground_roll_dist
from scipy.optimize import minimize, Bounds
import numpy as np
#E = 23652.9
#F = -22.93379
#Ga = 0.295879
#Gb = -5.46e-6
#H = 0
#
#N = 2 #engines
#
#h = 0 #ft
#T = 15 #C
#p = 1013.2 #HPa
#
#W = 133400 #lb Stage Length = 1
#
#C = .394884 # Flaps = 1+F
#
##Step 1
#Vc = 0
#cnt = corrected_net_thrust.value(E, F, Vc, Ga, h, Gb, H, T) #lbf
#
##Step 2
#
#Vcto = CAS_TakeOff.value(C, W)
#cnt = corrected_net_thrust.value(E, F, Vcto, Ga, h, Gb, H, T) #lbf
#
#B8 = 0.007701 #ft/lb
#theta = (T + 273)/288.15
#delta = 1
#
#S_TO8 = TO_ground_roll_dist.value(B8, theta, W, delta, N, cnt)

#Optimization

E = 23652.9
F = -22.93379
Ga = 0.295879
Gb = -5.46e-6
H = 0

N = 2 #engines

h = -11 #ft
T = 5.4 #C
p = 1010.0 #HPa
p0 = 1013.25 

W = 133400 #lb Stage Length = 1

C = .394884 # Flaps = 1+F

B8 = 0.007701 #ft/lb
theta = (T + 273)/288.15
delta = p/p0

d_r = 2535.904

cost = opt_TO_ground_roll_dist.cost

x0 = [139200, .95]
bounds = Bounds([.5*x0[0], .5], [1.5*x0[0], 1])
res = minimize(cost, x0, (B8, theta, delta, E, F, C, Ga, Gb, h, T, H, d_r),
        bounds=bounds)
