import numpy as np
from modules.ANP_steps import sin_of_climb_angle
import matplotlib.pyplot as plt
weights = np.array([133400, 139200, 145200, 155900, 169800])

N = 2 
Fn_delta = 20871
weights = weights / 0.9
R = 0.071

Vc = 150 

sins_gamma = sin_of_climb_angle(N, Fn_delta, weights, R, Vc)

plt.plot(weights, sins_gamma, 'o')
plt.xlabel('weight')
plt.ylabel('sin_gamma')
plt.show()
