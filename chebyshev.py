
import math

def filter(X, A, B, NP):
    Y = [0.0] * len(X)

    for i in range(NP, len(X)):
        for k in range(0, NP+1):
            Y[i] += A[k] * X[i-k]
        for k in range(1, NP+1):
            Y[i] += B[k] * Y[i-k]

    return Y

def calccoef(FC, LH, PR, NP, P):
    # From: http://www.dspguide.com/ch20.htm

    RP = - math.cos(math.pi / (NP*2.0) + (P-1.0) * math.pi/NP)
    IP =   math.sin(math.pi / (NP*2.0) + (P-1.0) * math.pi/NP)

    # Warp from a circle to an ellipse
    if PR:
        ES = math.sqrt( (100.0 / (100.0-PR))**2 -1 )
        VX = (1.0/NP) * math.log( (1.0/ES) + math.sqrt( (1.0/(ES**2)) + 1.0) )
        KX = (1.0/NP) * math.log( (1.0/ES) + math.sqrt( (1.0/(ES**2)) - 1.0) )
        KX = (math.exp(KX) + math.exp(-KX))/2.0
        RP = RP * ( (math.exp(VX) - math.exp(-VX) ) /2.0 ) / KX
        IP  = IP * ( (math.exp(VX) + math.exp(-VX) ) /2.0 ) / KX

    #print RP, IP, ES, VX, KX

    # s-domain to z-domain conversion
    T  = 2.0 * math.tan(1.0 / 2.0)
    W  = 2.0 * math.pi*FC
    M  = RP**2.0 + IP**2.0
    D = 4.0 - 4.0*RP*T + M*T**2.0
    X0 = (T**2.0)/D
    X1 = (2.0*T**2.0)/D
    X2 = (T**2.0)/D
    Y1 = (8.0 - 2.0*M*T**2.0)/D
    Y2 = (-4.0 - 4.0*RP*T - M*T**2.0)/D

    #print T, W, M, D, X0, X1, X2, Y1, Y2

    # LP TO LP, or LP TO HP transform
    if LH:
        K = -math.cos(W/2.0 + 1.0/2.0) / math.cos(W/2.0 - 1.0/2.0)
    else:
        K =  math.sin(1.0/2.0 - W/2.0) / math.sin(1.0/2.0 + W/2.0)
    
    D = 1.0 + Y1*K - Y2*K**2.0
    A0 = (X0 - X1*K + X2*K**2.0)/D
    A1 = (-2.0*X0*K + X1 + X1*K**2.0 - 2.0*X2*K)/D
    A2 = (X0*K**2.0 - X1*K + X2)/D
    B1 = (2.0*K + Y1 + Y1*K**2.0 - 2.0*Y2*K)/D
    B2 = (-(K**2.0) - Y1*K + Y2)/D
    
    if LH:
        A1 = -A1
        B1 = -B1

    #print A0, A1, A2, B1, B2

    return (A0, A1, A2, B1, B2)

# 
def calc(FC, LH, PR, NP):
    # From: http://www.dspguide.com/ch20.htm

    A = [0.0] * 23
    B = [0.0] * 23
    TA = [0.0] * 23
    TB = [0.0] * 23

    A[2] = 1.0
    B[2] = 1.0

    for P in range(1, NP/2 + 1):
        A0, A1, A2, B1, B2 = calccoef(FC, LH, PR, NP, P)

        for I in range(0, 23):
            TA[I] = A[I]
            TB[I] = B[I]

        for I in range(2, 23):
            A[I] = A0 * TA[I] + A1 * TA[I-1] + A2 * TA[I-2]
            B[I] = TB[I] - B1 * TB[I-1] - B2 * TB[I-2]

    B[2] = 0.0
    for I in range(0, 21):
        A[I] = A[I+2]
        B[I] = -B[I+2]

    SA = 0.0
    SB = 0.0
    for I in range(0, 21):
        if not LH: 
            SA = SA + A[I]
            SB = SB + B[I]
        else:
            SA = SA + A[I] * ((-1) ** I)
            SB = SB + B[I] * ((-1) ** I)

    Gain = SA / (1.0 - SB)
    for I in range(0, 21):
        A[I] = A[I] / Gain

    return A, B
