# RBR*argo*_DMQC
This repository contains the code to post-process data from RBR*argo*|2k CTDs on Argo floats. Users are invited to consult _Argo Quality Control Manual
For CTD and Trajectory Data_ (DOI: [10.13155/33951](http://dx.doi.org/10.13155/33951)), where more details are provided on how to properly use these tools.

These DMQC tools are provided "as is" to the community and should be used with caution and knowledge of their requirements. Both Python and MATLAB codes are made available.

### Compressibility_correction
For RBR*argo*|2k CTDs calibrated before April 2021, an additional correction must be applied to PSAL to account for compressibility errors. This repository includes:
* A database _RBRargo3_compressibility_table.csv_ listing the affected RBRargo|2k CTDs, along with the new calibration coefficients to be used. 
* A routine _RBRargo_compressibility_ that relies on the database mentioned above to automatically apply the compressibility correction if necessary.
* A routine _RBRargo_compressibility_example_ that shows an example on how _RBRargo_compressibility_ can be used.

### TEMP_CNDC
The variable TEMP_CNDC on RBR*argo*|2k CTDs corresponds to the temperature measured internally to the conductivity cell. It is a necessary variable to be able to correct for the long-term thermal inertia of the conductivity cell. If the variable is non-present, or the data quality is bad, TEMP_CNDC can be infered using the temperature history measure by the marine thermistor (TEMP variable) along with an IIR filter.
Repository includes:
* A routine _RBRargo3_TEMP_CNDC_from_TEMP_ADJUSTED_ that reconstructs the TEMP_CNDC using a recursive filter described in Lueck and Picklo (1990).
* A routine _RBRargo3_TEMP_CNDC_from_TEMP_ADJUSTED_example_ that shows an example on how _RBRargo3_TEMP_CNDC_from_TEMP_ADJUSTED_ can be used.

### Thermal inertia
Thermal inertia errors in conductivity measurements can be corrected using a combination of corrective algorithms, extensively described in Dever, Owens, Richards, Wijffels, Wong, Shkvorets, Halverson, and Johnson (2022)(DOI: [10.1175/JTECH-D-21-0186.1](http://dx.doi.org/10.1175/JTECH-D-21-0186.1))
Repository includes:
* A routine _RBRargo3_celltm_ that computes the thermal error associated with the thermal inertia of the conductivity cell (TEMP_celltm)
* A routine _RBRargo3_celltm_example_ that shows an example on how _RBRargo3_celltm_ can be used.

### SQRT_adjustment
In January 2024, RBR discovered an inconsistency between salinity readings reported by an instrument while streaming and the salinity readings calculated by Ruskin with the downloaded data ([RBR Report RBR#0016488revA](https://oem.rbr-global.com/floats/files/5898249/84673381/1/1710796993730/0016488revA+RBRargo3+C.T.D%2C+Field+Service+Bulletin+Q1+2024+bis.pdf)). The root cause of the discrepancy between the onboard calculated and the Ruskin calculated salinity was quickly attributed to an approximation in the square root calculation onboard the instrument. The error can be adjusted for in post-processing using the set of routines provided in the SQRT_adjustment folder. Repository includes:
* A routine _correct_sqrt_error_ that adjusts the salinity data for the square root discrepancy
* A routine _correct_sqrt_error_testcase_ that shows an example on how _correct_sqrt_error_ can be used.