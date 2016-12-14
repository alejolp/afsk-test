
import math
import random
import numpy as np

from matplotlib import pyplot as plt
from scipy.signal import butter, lfilter

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=2):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=2):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

SAMP = 9600 * 40
f1 = 9600
fcut = 9600 
brate = 9600

L = [-math.cos(2 * math.pi * f1 * (x) / SAMP) + math.cos(2 * math.pi * f1 * (x) / SAMP + math.pi) for x in range(2 * SAMP / brate)]

if False:
    Lx = [(1.0 * x / SAMP) for x in range(len(L))]
    plt.plot(L)
    plt.grid()
    plt.show()

A = [1.0 * math.cos(2 * math.pi * f1 * (x) / SAMP) for x in range(1 * SAMP / f1)]
As = [1.0 * math.sin(2 * math.pi * f1 * (x) / SAMP) for x in range(1 * SAMP / f1)]
B = [1.0 * math.cos(2 * math.pi * f1 * (x) / SAMP + math.pi) for x in range(1 * SAMP / f1)]

C = A * 4 + B + A + B * 8

#plt.plot(C, label='C')

dd = L.index(max(L))
#dd = SAMP / f1
#dd = int(0.0025 * SAMP)
#dd = int(round(0.000448 * SAMP))
#dd = int(round(0.00018229 * SAMP))
#dd = int(dd * 2.0 / 3.0)
#dd = int(SAMP * kk / 1e4)
#dd = kk
#dd = 31 # 8 20 31 70

D = [C[i] * C[(i + dd) % len(C)] for i in range(len(C))]
#D = [(C[i]**2) for i in range(len(C))]
#D = [(C[i]**2) * math.sin(2 * math.pi * (2*f1) * (i) / SAMP) for i in range(len(C))]
#D = [(C[i]**2) * math.cos(2 * math.pi * (2*f1) * (i) / SAMP) for i in range(len(C))]
#E = [(C[i]**2) * math.cos(2 * math.pi * (2*f1) * (i) / SAMP + math.pi) for i in range(len(C))]

if False:
    D = []
    for i in range(len(C)-len(A)):
        a = 0
        b = 0
        for k in range(len(A)):
            a += A[k] * C[k+i]
            b += As[k] * C[k+i]
        #D.append(a*a + b*b)
        D.append(math.atan2(b, a))

plt.plot(D, label='D')
#plt.plot(E, label='E')

#b, a = butter(6, (f1 / 2.0) / (0.5 * SAMP), btype='low', analog=False)
#E = lfilter(b, a, D)
#E = butter_bandpass_filter(D, fcut * 0.9, fcut * 1.1, SAMP)
E = butter_lowpass_filter(D, fcut, SAMP)

plt.plot(E, label="E, dd: {} ({})".format(dd, 1.0 * dd / SAMP))

# 
F = []
last = 0.5
th1 = 0.01
th2 = -0.01
for x in E:
    if x > th1:
        last = 0.5
    elif x < th2:
        last = -0.5
    F.append(last)

plt.plot(F, label='F')

plt.plot(A * (len(D) / len(A)), label='Clock')

plt.grid()
plt.legend()
plt.show()
