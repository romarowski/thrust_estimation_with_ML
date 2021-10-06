import numpy as np
import pdb
def sin_of_climb_angle(N, Fn_delta, W_delta, R, Vc, e=0 ):
    #Equation B-12 sin(gamma)
    #K: speed dependent constant
    #N: nbr of engines
    #Fn_delta: Mean Thrust
    #W_delta: Weight corrected by pressure
    #R: Drag over lift- dependent on flap setting
    #e: Bank angle, (epsilon)
    #Vc: CAS
    #e = 0  for now
    K = .95 if Vc > 200 else 1.01 
    sin_gamma = K * (N * (Fn_delta / W_delta) - R / np.cos(e))
    return sin_gamma

def TO_ground_roll_distance(B8, theta, W_delta, N, Fn_delta):
    #Equation B-9
    
    return B8*theta*(W_delta)**2/(N*Fn_delta)

def first_climb_CAS(C, W):
    #Equation B-15
    return C * np.sqrt(W)

def distance_accelerate(V_T2, V_T1, N, Fn_delta, W_delta, R, e, ROC):
    #Equation B-17

    f = 0.95 #Factor to account for 8kt headwind
    g = 32.17 #ft/s Gravity accel 
    k = 1.688 #ft/s per kt  Convert kt to ft/s
    e = 0 #Bank angle set to zero for now
    min2sec = 60 #sec/min
    V_T = (V_T2 + V_T1) / 2

    a_max = g * (N * (Fn_delta / W_delta) - R / np.cos(e))
    G = ROC / (min2sec * k * V_T) #Climb gradient

    s = f * k**2 * (V_T2**2 - V_T1**2) / 2 * (a_max - G * g)
                                            
    return s

def corrected_net_thrust(df_coeffs, Vc, T, h):
    E  = df_coeffs['E'].iloc[0]
    F  = df_coeffs['F'].iloc[0]
    Ga = df_coeffs['Ga'].iloc[0]
    Gb = df_coeffs['Gb'].iloc[0]
    H  = df_coeffs['H'].iloc[0] 
    return E + F*Vc + Ga*h + Gb*h**2 + H*T
                                               
