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

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def main():

    if "audio01.wav" in sys.argv[1]:
        pos_start = 267569
        pos_end = 309673
    elif "audio02.wav" in sys.argv[1]:
        pos_start = 325400
        pos_end = 363900

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

    E = butter_lowpass_filter(B, fcut, samp_rate, 3)

    # Schmitt trigger/hysteresis

    F = []
    if False:
        th1 = 65000
        th2 = -65000
    else:
        th1 = 4e6
        th2 = -4e6
    vh = th1
    vl = th2
    last = vh
    for x in E:
        if x > th1:
            last = vh
        elif x < th2:
            last = vl
        F.append(last)

    #plt.plot(A, label='A')
    plt.plot(B, label='B')
    plt.plot(E, label='E')
    plt.plot(F, label='F')
    #plt.plot([vh * math.sin(2 * math.pi * 1200 * x / samp_rate) for x in range(len(F))], label='Clock')
    plt.plot([(vh * (int((x + 10) * 1200 / samp_rate) % 2)) for x in range(len(F))], label='Clock', color='m')

    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
