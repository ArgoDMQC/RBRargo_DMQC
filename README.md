# README #

The repository includes resources necessary to improve data quality from the RBRargo3 in Delayed-Mode Quality Check (DMQC). It addresses two major aspects of the data:
1. Compressibility correction
2. Dynamic corrections

### Compressibility correction ###

Both MATLAB and python routines are provided to apply the compressibility corrections for the RBRargo3 CTDs deployed on Argo floats. Both routines are equivalent and complete the following steps:

1. Compute conductivity using gsw_C_from_SP( ) (GSW TEOS-10)
2. Update conductivity's compressibility calibration using

  *Cnew = compressibility_RBRargo3(COND,PRES,PRES_ADJUSTED,WMO)*

The routine automatically imports the updated calibration coefficients for the WMO numbers specified using the CSV file provided.

An example code is provided as a tarting point. For python users, a notebook is provided (*compressibility_correction_example.ipynb*), and a identical routine is provided for MATLAB users (*compressibility_correction_example.m*)



### Dynamic corrections ###

A single routine call *celltm_RBRargo3.m* does it all and apply all state-of-the-art corrections for the RBRargo3 CTD

* Requirements:
  * TEMP - Marine temperature measured by the RBRargo3 [in °C]
  * PRES - Sea Pressure measured by the RBRargo3 [in dbar]
  * TEMP_CNDC - Internal temperature recorded by the RBRargo3 [in °C]
  * e_time - the elapsed time [in s]

The input vectors should be organized in ascending order, that is where e_time is increasing.

1. Compute conductivity using gsw_C_from_SP( )
2. Compute elapsed time. If not directly available, it can be inferred from the sampling frequency (and NB_SAMPLE_CTD if data is binned).
3. Apply dynamic correction:

 Run *TEMPcell = celltm_RBRargo3(TEMP_ADJUSTED ,PRES_ADJUSTED,TEMP_CDNC, e_time)*
3. Compute corrected salinity using gsw_SP_from_C( ) with TEMPcell


### Who do I talk to? ###
Write to Mathieu Dever at argo@rbr-global.com
