# https://www.youtube.com/watch?v=JUqrVL6lOZ8 matlab/octave lesson

# Purpose: undertanding Fourier transform

import numpy as np
import matplotlib.pyplot as plt

# angular frequency of a chosen frequency
freq1 = 50.
w1 = 2.*np.pi*freq1

# angular frequency of a chosen frequency
freq2 = 5.
w2 = 2.*np.pi*freq2

# angular frequency of a chosen frequency
freq3 = 150.
w3 = 2.*np.pi*freq3

# set maximum time (in seconds)
maxt = 1.

# set sampling frequency
samplingfrequency = 1.e3 # in Hz
samples = int(maxt*samplingfrequency)
print('number of samples',samples)

# create time-axis
t = np.linspace(0.,maxt,samples)

# create periodic signal
signal1 = np.sin(w1*t)
signal2 = np.sin(w2*t)
signal3 = np.sin(w3*t)

# combine signals
signalfinal = signal1+signal2+signal3

# plot results and Fourier transform

fig=plt.figure(1,figsize=(8,8))
#plt.axis([xmin,xmax,ymin,ymax])
plt.subplot(211)
plt.plot(t,signalfinal)
plt.xlabel('time [s]')
plt.ylabel('amplitude')

a = np.fft.fft(signalfinal)
xf = np.linspace(0.,samplingfrequency,samples)
plt.subplot(212)
nyquist = samplingfrequency/2.
plt.axis([0.,nyquist,0.,np.max(a)])
plt.plot(xf,a)
plt.xlabel('frequency [Hz] max is Nyquist - half of sampling freq')
plt.ylabel('FFT power spectrum')

plt.show()
plt.close()