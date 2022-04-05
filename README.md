# Statistical Thrust Estimation Using ADS-B Data, Machine Learning, Optimization, Aerodynamics and the Aircraft Noise Performance Database

*Abstract:* Current aircrat operations noise modelling relies on the use of pre-defined vertical profiles (altitude vs. distance) to produce noise contours.
Altough computationally inexpensive this approach leads to noise contours that are not based in actual flight trajectories.
In this work we propose the use of real vertical profiles (obtained from ADS-B data) to obtain aircraft noise modelling parameters (like thrust, weight, flap
configuration, etc.) by performing an advanced curve fitting algorigthm similar to that introduced by Yang et. al. on Mathematical Programming for Piecewise 
Linear Regression Analysis (2016). We run a multi-objective optimization on the aerodynamic equations presented in Appendix B of ECAC, Doc 29: 4th Edition, Volume 2, 2016.
where the loss is the quadratic difference between the estimated vertical profile and Calibrated Airspeed vs the ADS-B measured profile and Calibrated airspeed. 
Two key contributions of this work are:
1. Robust identifation of aircraft *step type* (aircraft is climbing, accelerating, levelling or descending) solely from noisy ADS-B data.
2. The NoiseCraft.py library, that allows to clean and segment ADS-B data and estimate relevant paramters using high-level functions like:
        ```python:
        from modules import NoiseCraft
        craft = NoiseCraft(equipment='A320-211')
        craft.load_radar('radar.csv')
        craft.load_meteo('meteo.csv')
        craft.calculate_CAS()
        craft.segment()
        craft.recognize_steps()
        craft.new_vert_profile()
        ```


Algorithmic pipeline:
![pipeline](https://user-images.githubusercontent.com/36279027/161813955-475d2b74-09bb-4acb-80d0-ee7c673d1d46.png)

Some results:

![solution](https://user-images.githubusercontent.com/36279027/161814117-4ac1f07a-d449-4b3d-8fd3-6200073d502e.png)


This work was presented on the 35th AIRMOD Noise Specialist meeting on October 21’

**Keywords**: Machine Learning, Akaike Information Criterion, Piecewise Linear Regression, Python, Aviation, Noise Modelling, Pandas.
