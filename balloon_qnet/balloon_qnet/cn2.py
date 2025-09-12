import numpy as np

def hufnagel_valley(h, u_rms, Cn0):
    cn2 = 0.00594*((u_rms/27)**2)*((h*10**(-5))**10)*np.exp(-h/1000) + (2.7*10**(-16))*np.exp(-h/1500) + Cn0*np.exp(-h/100)
    return cn2

def slc_day(h):
    if np.size(h) == 1:
        h = np.ones(1)*h
    cn2 = np.zeros(np.size(h))
    layer_1 = h < 18.5
    layer_2 = h >= 18.5 and h < 240
    layer_3 = h >= 240 and h < 880
    layer_4 = h >= 880 and h < 7200
    layer_5 = h >= 7200 and h < 20000
    cn2[layer_1] = 1.7e-14
    cn2[layer_2] = 3.13e-13/h[layer_2]**1.05
    cn2[layer_3] = 1.3e-15
    cn2[layer_4] = 8.87e-7/h[layer_4]**3
    cn2[layer_5] = 2e-16/h[layer_5]**(1/2)
    return cn2

def slc_night(h): 
    cn2 = np.zeros(np.size(h))
    layer_1 = h < 18.5
    layer_2 = h >= 18.5 and h < 110
    layer_3 = h >= 110 and h < 1500
    layer_4 = h >= 1500 and h < 7200
    layer_5 = h >= 7200 and h < 20000
    cn2[layer_1] = 8.4e-15
    cn2[layer_2] = 2.87e-12/h[layer_2]**2
    cn2[layer_3] = 2.5e-16
    cn2[layer_4] = 8.87e-7/h[layer_4]**3
    cn2[layer_5] = 2e-16/h[layer_5]**(1/2)
    return cn2
