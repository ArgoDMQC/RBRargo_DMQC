#!/usr/bin/env python3
# coding: utf-8

def RBRargo3_TEMP_CNDC_from_TEMP_ADJUSTED(TEMP_ADJUSTED, elptime):
    """
    DESCRIPTION: This function infers the internal temperature TEMP_CNDC from
    the measured temperature TEMP_ADJUSTED profile using a Infinite Impulse
    response filter as defined in  Lueck and Picklo (1990), with alpha = 1.065,
    and tau = 154.8 s. TEMP_CNDC is used to apply the long-term thermal inertia
    correction to PSAL (see RBRargo3_celltm.m).

    It relies on an regularly sampled time series.


    I/O:
    Inputs:
    TEMP_ADJUSTED: temperature [°C] (ITS-90), as reported by the Argo float (size [mx1]) 
    elptime: elapsed time of the samples [seconds] (size [mx1])

    Outputs:
    TEMP_CNDC: Internal temperature [°C] infered from the TEMP_ADJUSTED measurements (size [mx1])


    !! Use this function at your own risk!!

    DEPENDENCIES:


    BIBLIOGRAPHY:
    Lueck, R. G., & Picklo, J. J. (1990). Thermal Inertia of Conductivity Cells:
    Observations with a Sea-Bird Cell, Journal of Atmospheric and Oceanic
    Technology, 7(5), 756-768

    Morison, J., Andersen, R., Larson, N., D'Asaro, E., & Boyd, T. (1994). The
    Correction for Thermal-Lag Effects in Sea-Bird CTD Data, Journal of
    Atmospheric and Oceanic Technology, 11(4), 1151-1164

    AUTHOR:
    Mathieu Dever, Qi Wang (e-mail: argo@rbr-global.com)
    20 Mar 2023

    v1.0 - 20/03/2023
    """

    import numpy as np
    import pandas as pd

    # Check inputs are vectors
    if TEMP_ADJUSTED.ndim != 1 or elptime.ndim != 1:
        raise Exception("Vector inputs only, please")

    # Check time is increasing
    if np.any(np.diff(elptime) < 0):
        raise Exception('Time is not increasing. Try sorting inputs first')
    elif np.any(np.diff(elptime) == 0) or len(np.unique(elptime)) != len(elptime):
        raise Exception('times are not unique')        

    # Set coefficients
    alpha = 1.065
    tau = 154.8
    # Sampling frequency
    fs = 1

    # Compute Tanomaly (long-term thermal mass temperature anomaly)

    # Interpolate needed variables onto a 1Hz time series
    gooddata = np.flatnonzero(~np.isnan(TEMP_ADJUSTED) & ~np.isnan(elptime))
    gridded = {}
    gridded["elptime"] = np.arange(min(elptime), max(elptime) + 2)
    gridded["TEMP"] = np.interp(gridded["elptime"], elptime[gooddata], TEMP_ADJUSTED[gooddata])
    
    # Nyquist frequency
    fn = fs / 2

    # Lueck and Picklo (1991)'s coefficients (a,b)
    a = 4 * fn * alpha * tau / (1 + 4 * fn * tau)
    b = 1 - 2 * a / alpha

    Tanomaly = np.zeros_like(gridded["TEMP"])
    for ii in range(1, len(gridded["TEMP"])):
        Tanomaly[ii] = -b * Tanomaly[ii - 1] + a * (gridded["TEMP"][ii] - gridded["TEMP"][ii - 1])

    TEMP_CNDC = gridded["TEMP"] - Tanomaly
    TEMP_CNDC = np.interp(elptime, gridded["elptime"], TEMP_CNDC)
    
    return TEMP_CNDC
