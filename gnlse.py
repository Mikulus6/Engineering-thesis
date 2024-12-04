# translated from matlab: https://github.com/jtravs/SCGBookCode/blob/master/gnlse.m

import numpy as np
from numpy.fft import *
from math import pi, log, factorial
from scipy.integrate import solve_ivp

eps = 2.2204e-16

last_printed_report_percentage = 0
def report(percentage):
    percentage = round(percentage * 100, 1)
    global last_printed_report_percentage
    if percentage > last_printed_report_percentage:
        last_printed_report_percentage = percentage
        print(f"{percentage} % complete")

def gnlse(T, A, w0, gamma, betas, loss, fr, RT, flength, nsaves):

    n = len(T)
    dT = T[2]-T[1]  # grid parameters
    V = 2*pi*np.arange(-n/2, n/2)/(n*dT)            # frequency grid
    alpha = log(10**(loss/10))                      # attenuation coefficient

    B = sum(betas[i]/factorial(i+2)*V**(i+2) for i in range(len(betas))) # Taylor expansion of betas

    L = 1j*B - alpha/2                              # linear operator

    if abs(w0) > eps :                              # if w0>0 then include shock
        gamma = gamma/w0
        W = V + w0                                  # for shock W is true freq
    else:
        W = 1                                       # set W to 1 when no shock

    RW = n*ifft(fftshift(RT))                       # frequency domain Raman
    L = fftshift(L)
    W = fftshift(W)                                 # shift to fft space

    def rhs(z, AW):
        report(z/flength)
        AT = fft(AW*np.exp(L*z))                    # time domain field
        IT = abs(AT)**2                             # time domain intensity
        if (len(RT) == 1) or (abs(fr) < eps):       # no Raman case
            M = ifft(AT*IT)                         # response function
        else:
            RS = dT*fr*fft(ifft(IT)*RW)             # Raman convolution
            M = ifft(AT*((1-fr)*IT + RS))           # response function
        R = 1j*gamma*W*M*np.exp(-L*z)               # full RHS of Eq. (3.13)
        return R

    Z = np.linspace(0, flength, nsaves)
    sol = solve_ivp(rhs, [0, flength], ifft(A), t_eval=Z, method='RK45')

    AW, AT = [],  []

    for i in range(nsaves):
        AW.append(sol.y[:, i] * np.exp(L * Z[i]))
        AT.append(fft(AW[i]))
        AW[i] = fftshift(AW[i] * dT * n)

    AT, AW = np.array(AT), np.array(AW)
    W = V + w0

    return Z, AT, AW, W
