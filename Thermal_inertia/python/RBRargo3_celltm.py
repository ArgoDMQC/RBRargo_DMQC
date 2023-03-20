#!/usr/bin/env python3
# coding: utf-8    

def RBRargo3_celltm(TEMP, PRES, TEMP_CNDC, e_time):
    
    """"
    DESCRIPTION: This function completes three thermal mass adjustements:

    * the thermal mass of the marine thermistor: a C-T lag is applied to
    re-align conductivity and temperature measurements and minimizes salinity
    spiking. Inputs are TEMP and e_time, output is Tcor (the aligned
    temperature data)

    * a long-term adjustement for thermal mass of the conductivity cell: this
    temperature anomaly (Tlong) is estimated using the internal temperature of
    the conductivity cell (TEMP_CNDC). Input is TEMP and TEMP_CNDC, output is
    Tlong.

    * a short-term adjustement for thermal mass of the conductivity cell: this
    temperature anomaly (Tshort) is estimated using the recursive algorithm
    developed by Lueck and Picklo (1990), and refined by Morison et al (1994).
    It requires two inputs: alpha, and tau.

    The appropriate coefficients for the correction are estimated using the
    profiling speed estimated from PRES and e_time.

    Data will be interpolated to 1 Hz for the correction.

    Note that the ucertainty in the correction is about the size of the
    correction itself.  This uncertainty should be treated as an independent
    error from other sources of error (such as the manufacturer's calibration
    uncertainty  and/or WJO/BS conductivity slope adjustment from historical
    T-S relationships) and thus should be combined in quadrature (square any
    relevant error terms, sum these squares, and then take the square root of
    the sum) for inclusion in psal_adjusted_error.

    I/O:
    Inputs: 
    TEMP: temperature [°C] (ITS-90), as reported by the Argo float (size [mx1]) 
    PRES: pressure [dbar] (size [mx1]) 
    TEMP_CNDC: Internal temperature [°C] reported by the RBRargo3 CTD (size [mx1])
    e_time: elapsed time of the samples [seconds] (size [mx1])

    Outputs: 
    TEMPcell: Temperature adjusted for thermal mass errors (size [mx1]).

    while TEMP is still the appropriate water temperature, TEMPcell should be used
    to adjust Salinity. (e.g. PSAL_ADJUSTED = gsw_SP_from_C(C0,TEMPcell,PRES)

    !! Use this function at your own risk!!
    The author takes no responsibility for errors or omissions, but will be
    happy to receive suggestions for improvements or corrections so that they
    can be reviewed, implimented if warranted, and passed on to others.  Users
    should e-mail the author so that they can receive notice of any updates.

    DEPENDENCIES:

    BIBLIOGRAPHY:
    Dever M., B. Owens, C. Richards, S. Wijffels, A. Wong, I. Shkvorets, M.
    Halverson, and G. Johnson (in prep.).﻿Static and dynamic performance of the
    RBRargo3

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
    if TEMP.ndim == 1:
        pass
    else:
        raise Exception("Vector inputs only, please")
    
    
    # Check time is increasing
    if any(np.diff(e_time) < 0):
        raise Exception('Time is not increasing. Try sorting inputs first')
    elif any(np.diff(e_time) == 0) or len(np.unique(e_time)) != len(e_time):
        raise Exception('times are not unique')
        
    # Compute the magitude of the ascent rate Vp
    Vp = np.zeros_like(TEMP)
    # Use centre-differencing
    Vp[1:-1] = np.abs((PRES[2:] - PRES[:-2]) / (e_time[2:] - e_time[:-2]))
    # use forward differencing for 1st datapoint to not lose it
    Vp[0] = np.abs((PRES[1] - PRES[0]) / (e_time[1] - e_time[0]))
    # use backward differencing for last datapoint to not lose it
    Vp[-1] = np.abs((PRES[-1] - PRES[-2]) / (e_time[-1] - e_time[-2]))
    # convert ascent rate to cm/s
    Vp *= 100
    
    # Compute thermal mass coefficients based on Vp.
    ctcoeff = 0.14 * Vp ** (-1.00)
    alpha = 0.37 * Vp ** (-1.03)
    tau = 16.02 * Vp ** (-0.26)
    
    # Compute Tcor (C-T lag)
    CTlag = -0.35
    Tcor = np.interp(e_time - CTlag, e_time, TEMP)
    
    # Compute Tlong (long-term thermal mass temperature anomaly)
    Tlong = ctcoeff * (TEMP_CNDC - Tcor)
    
    # Compute Tshort (long-term thermal mass temperature anomaly)
    # Interpolate needed variables onto a 1Hz time series
    gooddata = np.flatnonzero(~np.isnan(TEMP) & ~np.isnan(e_time))
    gridded = {}
    gridded["e_time"] = np.arange(np.min(e_time), np.max(e_time) + 2)
    gridded["TEMP"] = np.interp(gridded["e_time"], e_time[gooddata], TEMP[gooddata])
    
    # Sampling frequency
    fs = 1
    # Nyquist frequency
    fn = fs / 2
    
    # Lueck and Picklo (1991)'s coefficients (a,b)
    a = 4 * fn * alpha * tau / (1 + 4 * fn * tau)
    b = 1 - 2 * a / alpha
    gridded["a"] = np.interp(gridded["e_time"], e_time[gooddata], a[gooddata])
    gridded["b"] = np.interp(gridded["e_time"], e_time[gooddata], b[gooddata])
    
    # Apply Lueck and Picklo (1991) recursive filter
    gridded["Tshort"] = np.zeros_like(gridded["TEMP"])

    # First value = 0°C, so loop from second value
    for tt in range(1, len(gridded["TEMP"])):
        gridded["Tshort"][tt] = -gridded["b"][tt] * gridded["Tshort"][tt-1] + gridded["a"][tt] * (gridded["TEMP"][tt]-gridded["TEMP"][tt-1])

    # Go back and set the first temperature adjustment to the second, since the
    # float is probably rising as it takes its first samples
    gridded["Tshort"][0]=gridded["Tshort"][1]

    # Re-grid onto original time grid
    Tshort = np.nan * TEMP
    Tshort[gooddata] = np.interp(e_time[gooddata],gridded["e_time"],gridded["Tshort"])
    
    # Compute TEMPcell
    # Compute the adjusted temperature that should be used to compute
    # salinity corrected for thermal mass errors. It should ONLY be used to
    # calculate the adjusted salinity, it does not reflect the true temperature
    # of the water. Note the signs, as per the convention established in
    # Morison (1994)

    TEMPcell = Tcor + Tlong - Tshort

    return TEMPcell
    