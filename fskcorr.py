
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

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

SAMP = 115200
f1 = 1200
f2 = 2200
fcut = 500

L = [-math.cos(2 * math.pi * f1 * (x) / SAMP) + math.cos(2 * math.pi * f2 * (x) / SAMP) for x in range(2 * SAMP / 1200)]

if True:
    Lx = [(1.0 * x / SAMP) for x in range(len(L))]
    plt.plot(L)
    plt.grid()
    plt.show()

A = [1.0 * math.cos(2 * math.pi * f1 * (x) / SAMP) for x in range(10 * SAMP / f1)]
B = [1.0 * math.cos(2 * math.pi * f2 * (x) / SAMP) for x in range(20 * SAMP / f2)]

if False:
    C = [(x + 0.2 * (random.randint(0, 255) - 128) / 128.0) for x in (A + B + A + B)]
else:
    C = (A + B + A + B)

plt.plot(C)

dd = L.index(max(L))
#dd = int(0.0025 * SAMP)
#dd = int(round(0.000448 * SAMP))
#dd = int(round(0.00018229 * SAMP))
#dd = int(dd * 2.0 / 3.0)
#dd = int(SAMP * kk / 1e4)
#dd = kk
#dd = 31 # 8 20 31 70

D = [C[i] * C[(i + dd) % len(C)] for i in range(len(C))]

# plt.plot(D)

#b, a = butter(6, (f1 / 2.0) / (0.5 * SAMP), btype='low', analog=False)
#E = lfilter(b, a, D)
E = butter_lowpass_filter(D, fcut, SAMP)

plt.plot(E, label="dd: {} ({})".format(dd, 1.0 * dd / SAMP))

# 
F = []
last = 0.5
th1 = 0.079
th2 = -0.079
for x in E:
    if x > th1:
        last = 0.5
    elif x < th2:
        last = -0.5
    F.append(last)

plt.plot(F, label='F')

plt.grid()
plt.legend()
plt.show()
