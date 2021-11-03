import numpy as np
import matplotlib.pyplot as plt

def main():
    N = 2
    K = 0.95
   # R = np.array([0.061662, 0.096267, 0.067463, 0.101204, 0.11586, 0.057558, 
   #               0.066827, 0.071403, 0.056281])
    R = np.array([0.066827, 0.071403, 0.056281])

    W  = 139200
    Fn = 23255

    bank_angle = 0*np.pi/180

    sin_gamma = K * (N * Fn / W - R / np.cos(bank_angle))
    plt.subplot(211)
    plt.plot(R, sin_gamma, 'o')
    plt.xlabel('R')
    plt.ylabel('sin$\gamma$')

    dists = 1000 / np.tan(np.arcsin(sin_gamma))
    
    plt.subplot(212)
    plt.plot(R, dists, 'o')
    plt.xlabel('R')
    plt.ylabel('distances')
    plt.show()
    
    return 0 

if __name__ == '__main__':
    main()

