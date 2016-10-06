#!/usr/bin/env python
# coding: utf-8

# Schmitt trigger

import os, sys, math, wave, struct

import matplotlib.pyplot as plt
import numpy as np


def gen_sin(freq, samp_rate, samples, mult):
    return [mult * math.sin(2.0 * math.pi * freq * i / samp_rate) for i in range(samples)]

def gen_cos(freq, samp_rate, samples, mult):
    return [mult * math.cos(2.0 * math.pi * freq * i / samp_rate) for i in range(samples)]

def window_triangular(L):
    N = len(L)
    return [L[i] * (1.0 - abs((i - (N-1) / 2.0) / ((N-1)/2.0))) for i in range(len(L))]

def fmul(A, B, Astart = 0, Bstart = 0):
    R = 0
    for i in range(min(len(A), len(B))):
        R += A[(i + Astart) % len(A)] * B[(i + Bstart) % len(B)]
    return R

def avg(L):
    return sum(L) / len(L)

def median(L):
    M = sorted(L)
    return M[len(M) / 2]

def main():
    w = wave.open(sys.argv[1], 'rb')

    ch = w.getnchannels()
    samp_rate = w.getframerate()
    swidth = w.getsampwidth()

    print ch, samp_rate, swidth

    samples_per_bit = samp_rate / 1200

    if (samp_rate % 1200) > 0:
        samples_per_bit += 1

    if swidth == 2:
        fmtstr = "<h"
        sin_mult = 32767
    else:
        assert False

    print "samples_per_bit", samples_per_bit

    tas = gen_sin(1200, samp_rate, samples_per_bit, sin_mult) # * 0.83
    tac = gen_cos(1200, samp_rate, samples_per_bit, sin_mult)

    tbs = gen_sin(2200, samp_rate, samples_per_bit, sin_mult) # * 1.33
    tbc = gen_cos(2200, samp_rate, samples_per_bit, sin_mult)

    if False:
        tas = window_triangular(tas)
        tac = window_triangular(tac)
        tbs = window_triangular(tbs)
        tbc = window_triangular(tbc)

    if False:
        window = np.hanning(samples_per_bit)

        tas = tas * window
        tac = tac * window
        tbs = tbs * window
        tbc = tbc * window

    print tas

    samples_history = [([0] * samples_per_bit) for i in range(ch)]
    samples_history_idx = 0

    samples_result_a = [[0] for i in range(ch)]
    samples_result_b = [[0] for i in range(ch)]
    samples_result = [[0] for i in range(ch)]
    schmitt_result = [[0] for i in range(ch)]

    thresh_high = (sin_mult ** 2) / 20
    thresh_low = -(sin_mult ** 2) / 20

    #thresh_high = -8 + 14
    #thresh_low = -8 - 14

    end = False
    while not end:
        for chnum in range(ch):
            datastr = w.readframes(1)
            if len(datastr) < 1:
                end = True
                break
            data = struct.unpack(fmtstr, datastr)

            samples_history[chnum][samples_history_idx] = data[0]
            samples_history_idx = (samples_history_idx + 1) % len(samples_history[chnum])

            #this_slice = samples_history[chnum] 
            #this_slice = [samples_history[chnum][(i + samples_history_idx) % len(window)] for i in range(len(window))] * window

            #mula = math.sqrt(fmul(this_slice, tas) ** 2 + fmul(this_slice, tac) ** 2)
            #mulb = math.sqrt(fmul(this_slice, tbs) ** 2 + fmul(this_slice, tbc) ** 2)

            mula = math.sqrt(fmul(samples_history[chnum] , tas, samples_history_idx) ** 2 + fmul(samples_history[chnum] , tac, samples_history_idx) ** 2)
            mulb = math.sqrt(fmul(samples_history[chnum] , tbs, samples_history_idx) ** 2 + fmul(samples_history[chnum] , tbc, samples_history_idx) ** 2)

            #mula = 20 * math.log(mula, 10) if mula > 0 else 0
            #mulb = 20 * math.log(mulb, 10) if mulb > 0 else 0

            maxa = max(samples_result_a[chnum][-(samples_per_bit * 8):])
            maxb = max(samples_result_b[chnum][-(samples_per_bit * 8):])

            maxmid = (maxa + maxb) / 2

            agca = (maxmid / maxa) if (abs(maxa) > 0) else 1
            agcb = (maxmid / maxb) if (abs(maxb) > 0) else 1

            #mula *= agca
            #mulb *= agcb

            diff = (mula - mulb)

            thresh_high = maxmid / 4
            thresh_low = maxmid / -4

            if diff >= thresh_high:
                schmitt_value = thresh_high
            elif diff <= thresh_low:
                schmitt_value = 0
            else:
                schmitt_value = schmitt_result[chnum][-1]

            samples_result[chnum].append(diff)
            samples_result_a[chnum].append(mula)
            samples_result_b[chnum].append(mulb)
            schmitt_result[chnum].append(schmitt_value)

    print "len samples_result", len(samples_result[0]), max(samples_result[0]), min(samples_result[0])

    if True:
        plt.plot(samples_result[0], label="a - b")
        plt.plot(samples_result_a[0], label="a")
        plt.plot(samples_result_b[0], label="b")
        plt.plot(schmitt_result[0], label="sch")

#    plt.plot([max(samples_result_a[0][i:i+(samples_per_bit * 8)]) for i in range(len(samples_result_a[0]) - samples_per_bit * 8)], label="max(a)")
#    plt.plot([max(samples_result_b[0][i:i+(samples_per_bit * 8)]) for i in range(len(samples_result_b[0]) - samples_per_bit * 8)], label="max(b)")

    if False:
        plt.plot(tas)
        plt.plot(tac)
        plt.plot(tbs)
        plt.plot(tbc)

    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
