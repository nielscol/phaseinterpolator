from __future__ import division
from math import sqrt
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import numpy as np
from copy import copy


# Golden-section search 
def gss_find_arg(func, arg, target, params, _min, _max, conv_tol):
    params_a = copy(params)
    params_b = copy(params)
    phi = (1.0 + sqrt(5.0))/2.0
    iter = 1
    err = np.inf # Initial error (%)
    # Print column headers
    #delta = (_max-_min)*((phi-1)/(3-2*phi))
    #guess_h = _max + delta
    #guess_l = _min - delta
    guess_l = a = _min
    guess_h = b = _max
    delta = 0
    #print "#", delta, guess_h, guess_l
    # Iterate until termination criterion is reached
    while err > conv_tol:
        delta = (phi - 1.0)*(b - a)
        guess_a = a + delta
        guess_b = b - delta
        params_a[arg] = guess_a
        params_b[arg] = guess_b
        # print "%", guess_a, guess_b, abs(target - func(**params_a))**2, abs(target - func(**params_b))**2
        if abs(target - func(**params_b)) < abs(target - func(**params_a)):
            b = guess_a
            val = b
        else:
            a = guess_b
            val = a
        if val == 0.0:
            err = np.inf
        else:
            err = ((2.0 - phi)*abs((guess_a - guess_b)/val))
        params[arg] = val
        iter += 1
    return val

