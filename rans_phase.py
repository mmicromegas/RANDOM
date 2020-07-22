import numpy as np
import matplotlib.pyplot as plt
import random

# create periodic signal
freq = 5.
w = 2.*np.pi*freq

# set maximum length of domain
maxd = 1.

# set resolution "grid" points
sampling = 1.e3 
samples = maxd*sampling

# create time-axis
d = np.linspace(0.,maxd,samples)

ampl = 3.

# create periodic signal
signal = ampl*np.sin(w*d)

plt.figure(figsize=(7,6))
plt.plot(d,signal)

signal_shift = []

ns = 260
for i in range(ns):
    # randomly vary signals phase
    shift = random.random()*2.*np.pi
#    print(i,shift)	
    signal_shift.append(ampl*np.sin(w*d+shift)) 
    plt.plot(d,signal_shift[i])

signal_final = signal_shift[0] 	
for i in range(1,ns):
#    print(i)
    signal_final = signal_final*signal_shift[i]
	
plt.figure(figsize=(7,6))
plt.plot(d,signal_final/max(np.abs(signal_final)))	
	
plt.show(block=False)
