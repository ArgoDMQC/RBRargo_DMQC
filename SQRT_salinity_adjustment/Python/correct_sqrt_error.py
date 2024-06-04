#!/usr/bin/env python3
# coding: utf-8

def correct_sqrt_error(Sraw):
    """
    DESCRIPTION: This function post-corrects salinity measurements to account for a
    discontinuity in onboard salinity computation


    This function post-corrects salinity measurements to account for a
    discontinuity in onboard salinity computation
    
    VARIABLES:
    Sraw - measured salinity
    Sadj - adjusted salinity

    AUTHOR:
    Mathieu Dever, Qi Wang (e-mail: argo@rbr-global.com)
    09 May 2024

    v1.0 - 09/05/2024
    """

    import numpy as np

    # Adjusts salinity
    Sadj = np.copy(Sraw)

    # Set correction flag
    correction = False

    # Loops through the time series
    for idx, val in enumerate(Sraw):
        if correction == False:
            if val < 35.000:
                correction = True
        else:
            if val > 35.002:
                correction = False
        if correction == True:
            error = (3.559e-10)*np.exp(0.4403*val)
            Sadj[idx] = val - error
    
    return Sadj
