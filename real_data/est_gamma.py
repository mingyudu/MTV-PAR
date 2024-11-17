# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:26:26 2019

"""


import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd
#import math
#import sys
import scipy
from scipy import signal
#reference: https://papers.nips.cc/paper/6505-fast-active-set-methods-for-online-spike-inference-from-calcium-imaging.pdf
#github: https://github.com/j-friedrich/OASIS/blob/master/oasis/functions.py

def estimate_parameters(y, p=2, range_ff=[0.25, 0.5], method='mean', lags=10, fudge_factor=1., nonlinear_fit=False):
    """
    Estimate noise standard deviation and AR coefficients
    Parameters
    ----------
    p : positive integer
        order of AR system
    lags : positive integer
        number of additional lags where he autocovariance is computed
    range_ff : (1,2) array, nonnegative, max value <= 0.5
        range of frequency (x Nyquist rate) over which the spectrum is averaged
    method : string, optional, default 'mean'
        method of averaging: Mean, median, exponentiated mean of logvalues
    fudge_factor : float (0< fudge_factor <= 1)
        shrinkage factor to reduce bias
    """

    sn = GetSn(y, range_ff, method)
    g = estimate_time_constant(y, p, sn, lags, fudge_factor, nonlinear_fit)

    return g, sn


def estimate_time_constant(y, p=2, sn=None, lags=10, fudge_factor=1., nonlinear_fit=False):
    """
    Estimate AR model parameters through the autocovariance function
    Parameters
    ----------
    y : array, shape (T,)
        One dimensional array containing the fluorescence intensities with
        one entry per time-bin.
    p : positive integer
        order of AR system
    sn : float
        sn standard deviation, estimated if not provided.
    lags : positive integer
        number of additional lags where he autocovariance is computed
    fudge_factor : float (0< fudge_factor <= 1)
        shrinkage factor to reduce bias
    Returns
    -------
    g : estimated coefficients of the AR process
    """

    if sn is None:
        sn = GetSn(y)

    lags += p
    # xc = axcov(y, lags)[lags:]
    y = y - y.mean()
    xc = np.array([y[i:].dot(y[:-i if i else None]) for i in range(1 + int(lags))]) / len(y)

    if nonlinear_fit and p <= 2:
        xc[0] -= sn**2
        g1 = xc[:-1].dot(xc[1:]) / xc[:-1].dot(xc[:-1])
        if p == 1:
            def func(x, a, g):
                return a * g**x
            popt, pcov = curve_fit(func, list(range(len(xc))), xc, (xc[0], g1)) #, bounds=(0, [3 * xc[0], 1]))
            return popt[1:2] * fudge_factor
        elif p == 2:
            def func(x, a, d, r):
                return a * (d**(x + 1) - r**(x + 1) / (1 - r**2) * (1 - d**2))
            popt, pcov = curve_fit(func, list(range(len(xc))), xc, (xc[0], g1, .1))
            d, r = popt[1:]
            d *= fudge_factor
            return np.array([d + r, -d * r])

    xc = xc[:, np.newaxis]
    A = scipy.linalg.toeplitz(xc[np.arange(int(lags))],
                              xc[np.arange(int(p))]) - sn**2 * np.eye(int(lags), int(p))
    g = np.linalg.lstsq(A, xc[1:])[0]
    gr = np.roots(np.concatenate([np.array([1]), -g.flatten()]))
    gr = (gr + gr.conjugate()) / 2.
    gr[gr > 1] = 0.95 + np.random.normal(0, 0.01, np.sum(gr > 1))
    gr[gr < 0] = 0.15 + np.random.normal(0, 0.01, np.sum(gr < 0))
    g = np.poly(fudge_factor * gr)
    g = -g[1:]

    return g.flatten()

def GetSn(y, range_ff=[0.25, 0.5], method='mean'):
    """
    Estimate noise power through the power spectral density over the range of large frequencies
    Parameters
    ----------
    y : array, shape (T,)
        One dimensional array containing the fluorescence intensities with
        one entry per time-bin.
    range_ff : (1,2) array, nonnegative, max value <= 0.5
        range of frequency (x Nyquist rate) over which the spectrum is averaged
    method : string, optional, default 'mean'
        method of averaging: Mean, median, exponentiated mean of logvalues
    Returns
    -------
    sn : noise standard deviation
    """

    ff, Pxx = scipy.signal.welch(y)
    ind1 = ff > range_ff[0]
    ind2 = ff < range_ff[1]
    ind = np.logical_and(ind1, ind2)
    Pxx_ind = Pxx[ind]
    sn = {
        'mean': lambda Pxx_ind: np.sqrt(np.mean(Pxx_ind / 2)),
        'median': lambda Pxx_ind: np.sqrt(np.median(Pxx_ind / 2)),
        'logmexp': lambda Pxx_ind: np.sqrt(np.exp(np.mean(np.log(Pxx_ind / 2))))
    }[method](Pxx_ind)

    return sn

#y = np.genfromtxt("c:/study/y_spike.txt", delimiter=None)
#y = y[:,1]

#rst = estimate_parameters(y, p=1)






