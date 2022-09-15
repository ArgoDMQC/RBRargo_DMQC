#!/usr/bin/env python
# coding: utf-8

# In[71]:


def compressibility_RBRargo3(COND,PRES,PRES_ADJUSTED,WMO):
  
    """
    DESCRIPTION: This function updates the compressibility calibration on a set
    of Argo float quipped with RBRargo3 CTDs. It imports the new set of
    coefficients from a look-up table accompanying this routine.

    I/O:
    Inputs:
    COND: Conductivity [mS/cm], as reported by the Argo float (size [mx1])
    PRES: pressure [dbar] (size [mx1])
    WMO: WMO identifier of the float, required to find it in the look-up table

    Outputs:
    Cnew: Conductivity [mS/cm] updated with the newest compressibility
    calibration for the specific WMO.

    !! Use this function at your own risk!!

    DEPENDENCIES:
    RBRargo3_compressibility_table.csv

    BIBLIOGRAPHY:


    AUTHOR:
    Mathieu Dever (e-mail: argo@rbr-global.com)
    15 Sep 2022

    v1.0 - 15/09/2022
    """
    import numpy as np
    import pandas as pd

    # Checks
    # Check inputs are vectors
    if len(COND.shape)>1:
        raise Exception("Vector inputs only, please")

    # Check sizes are right
    if (COND.shape != PRES.shape) | (COND.shape != PRES_ADJUSTED.shape) | (PRES_ADJUSTED.shape != PRES.shape):
        raise Exception("Input size not matching")
    
    if np.size(WMO)!=1:
        raise Exception("WMO size must be 1")

    # Import look-up table with new coefficients      
    coefs = pd.read_csv('RBRargo3_compressibility_table.csv',
                        usecols=range(0,13))

    if len(coefs.X2[coefs.WMO == WMO])==0:
        raise Exception("WMO not found in the database")
    
    if any(np.isnan(coefs.X2[coefs.WMO == WMO])):
        raise Exception("Sorry, No coefficients found for this WMO.")
    
    
    X2old = np.array(coefs.X2old[coefs.WMO == WMO])
    X3old = np.array(coefs.X3old[coefs.WMO == WMO])
    X4old = np.array(coefs.X4old[coefs.WMO == WMO])
    X2 = np.array(coefs.X2[coefs.WMO == WMO])
    X3 = np.array(coefs.X3[coefs.WMO == WMO])
    X4 = np.array(coefs.X4[coefs.WMO == WMO])

    Cnew = COND*(1 + X2old*PRES + X3old*PRES**2 + X4old*PRES**3)/(1 + X2*PRES_ADJUSTED + X3*PRES_ADJUSTED**2 + X4*PRES_ADJUSTED**3)
    
    return Cnew

