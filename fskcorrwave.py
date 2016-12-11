#!/usr/bin/env python
# coding: utf-8

import math
import random
import wave
import struct
import sys
import os

import numpy as np

from matplotlib import pyplot as plt
from scipy.signal import butter, lfilter

import chebyshev

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_lowpass_filter2(data, cutoff, fs, order=2):
    A, B = chebyshev.calc(cutoff / (0.5 * fs), 0, 0.0, order)
    return chebyshev.filter(data, A, B, order)

def main():

    if "audio01.wav" in sys.argv[1]:
        pos_start = 267569
        pos_end = 309673
        wth1 = 1e7
        wth2 = -1e7
        th1 = 1e7
        th2 = -1e7
        phcorr = 0
    elif "audio02.wav" in sys.argv[1]:
        pos_start = 325400
        pos_end = 363900
        wth1 = 100000
        wth2 = -250000
        th1 = 6e6
        th2 = -3.5e6
        phcorr = -20

    w = wave.open(sys.argv[1], 'rb')

    ch = w.getnchannels()
    samp_rate = w.getframerate()
    swidth = w.getsampwidth()

    print ch, samp_rate, swidth

    samples_per_bit = int(math.ceil(samp_rate / 1200.0))

    if swidth == 2:
        fmtstr = "<" + ("h" * ch)
    else:
        assert False

    print "samples_per_bit", samples_per_bit

    bps_rate = 1200
    f1 = 1200
    f2 = 2200
    fcut = 600 # Hz for LPF

    L = [-math.cos(2 * math.pi * f1 * (x) / samp_rate) + math.cos(2 * math.pi * f2 * (x) / samp_rate) for x in range(2 * samp_rate / bps_rate)]
    dd = L.index(max(L))

    if False:
        Lx = [(1.0 * x / samp_rate) for x in range(len(L))]
        plt.plot(L)
        plt.grid()
        plt.show()

    A = []

    w.setpos(pos_start)

    end = False
    while not end and w.tell() < (pos_end):
        datastr = w.readframes(1)

        if len(datastr) < 1:
            end = True
            break

        data = struct.unpack(fmtstr, datastr)

        A.append(data[0])

    B = [A[i] * A[(i + dd) % len(A)] for i in range(len(A))]

    E = butter_lowpass_filter2(B, fcut, samp_rate)

    # Schmitt trigger/hysteresis

    if True:
        F = []
        vh = th1
        vl = th2
        last = vh
        for x in E:
            if x > th1:
                last = vh
            elif x < th2:
                last = vl
            F.append(last)
    elif False:
        hyst_hist = [0] * (samples_per_bit * 16)
        hyst_hist_idx = 0
        F = []
        vh = 0.0
        vl = 0.0
        last = 0
        for x in E:
            hyst_hist[hyst_hist_idx % len(hyst_hist)] = x
            hyst_hist_idx += 1
            L = sorted(hyst_hist)
            th1 = vh = L[9 * len(L) / 10] * 0.75
            th2 = vl = L[1 * len(L) / 10] * 0.75
            if x > th1:
                last = 1
            elif x < th2:
                last = 0
            F.append(last * wth1)
    else:
        F = []
        G = []

        eta = 1.6

        eta_fast = eta / 10.0
        mean_fast = 0.0
        median_fast = 0.0

        eta_slow = eta / 100.0
        mean_slow = 0.0
        median_slow = 0.0

        last = 0
        for x in E:
            mean_fast += eta_fast * (x - mean_fast)
            #median_fast += eta_fast * math.copysign(1.0, x - median_fast)

            mean_slow += eta_slow * (x - mean_slow)
            #median_slow += eta_slow * math.copysign(1.0, x - median_slow)

            d = mean_fast - mean_slow
            F.append(d)

            if d >= wth1:
                last = wth1
            elif d < wth2:
                last = wth2

            G.append(last)

    #plt.plot(A, label='A')
    
    #plt.plot(B, label='B')
    plt.plot(E, label='E')
    #plt.hist(E, bins=50, log=True, label='E_hist')
    #plt.plot(F, label='F')
    #plt.hist(F[1200:18500], bins=50, log=True, label='F_hist')
    #plt.hist(G, bins=50, label='G_hist')
    #plt.plot(G, label='G')
    #plt.plot([vh * math.sin(2 * math.pi * 1200 * x / samp_rate) for x in range(len(F))], label='Clock')
    
    plt.plot([(wth1 * (int((x + phcorr) * 1200 / samp_rate) % 2)) for x in range(len(F))], label='Clock', color='m')

    #plt.ylim(-5, 5)
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
