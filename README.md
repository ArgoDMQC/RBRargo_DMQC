# README #

The repository includes resources necessary to improve data quality from the RBRargo3 in Delayed-Mode Quality Check (DMQC). It addresses two major aspects of the data:
1. Compressibility correction for RBRargo3 CTDs calibrated before April 2021. For RBRargo3 CTDs calibrated after April 2021, the compressibility of the conductivity sensor is corrected for as part of the routine CTD calibration.
2. Dynamic corrections

### Compressibility correction ###

A single routine call *compressibility_RBRargo3.m* does it all and apply the compressibility corrections for the RBRargo3 CTDs deployed on Argo floats.

1. Compute conductivity using gsw_C_from_SP( )
2. Update conductivity's compressiblity calibration using

  *Cnew = compressibility_RBRargo3(COND,PRES,PRES_ADJUSTED,WMO)*

  The routine automatically imports the updated calibration coefficients for the WMO numbers specified.

3. Compute corrected salinity using gsw_SP_from_C( ) with Cnew


### Dynamic corrections ###

A single routine call *celltm_RBRargo3.m* does it all and apply all state-of-the-art corrections for the RBRargo3 CTD

* Requirements:
  * TEMP - Marine temperature measured by the RBRargo3 [in °C]
  * PRES - Sea Pressure measured by the RBRargo3 [in dbar]
  * NB_SAMPLE_CTD
  * TEMP_CNDC - Internal temperature recorded by the RBRargo3 [in °C]
  * sampling frequency [in Hz]

1. Compute conductivity using gsw_C_from_SP( )
2. Compute elapsed time. If not directly available, it can be inferred from the sampling frequency (and NB_SAMPLE_CTD if data is binned).
3. Apply dynamic correction:

 Run *TEMPcell = celltm_RBRargo3(TEMP_ADJUSTED ,PRES_ADJUSTED,TEMP_CDNC, elptime)*

4. Compute corrected salinity using gsw_SP_from_C( ) with TEMPcell


### Who do I talk to? ###
Write to RBR at argo@rbr-global.com
