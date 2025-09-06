function TEMPcell = RBRargo3_celltm(TEMP,PRES,TEMP_CNDC,e_time, RBRargo_CTcell_model)

%{
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
Mathieu Dever, Hanyuan Liu (e-mail: argo@rbr-global.com)
05 Sept 2025

v1.0 - 22/11/2021
v2.0 - 07/07/2022 - update the relationships between key coefficients and
water-speed (based on revisions to Dever et al. 2022)
v3.0 - 05/09/2025 - add option for RBRargo_CTcell_model
%}

%% Set default if argument not provided
if nargin < 5 || isempty(RBRargo_CTcell_model)
    warning('No RBRargo_CTcell_model argument provided, defaulting to "RBRargo|2k_CTcell_pre2025" coefficients.');
    RBRargo_CTcell_model = 'RBRargo|2k_CTcell_pre2025';
end

%% checks 

% Check inputs are vectors
[n,m]=size(TEMP);
if min([n,m])>1
    disp('Vector inputs only, please')
    return
end
clear m n

% Check time is increasing
if any(diff(e_time)<0)
    error('Time is not increasing. Try sorting inputs first')
elseif any(diff(e_time)==0) || length(unique(e_time))~=length(e_time)
    error('times are not unique')
end

%% Compute the magitude of the ascent rate Vp

Vp = zeros(size(TEMP));
% Use centre-differencing
Vp(2:end-1) = abs((PRES(3:end)-PRES(1:end-2))./(e_time(3:end)-e_time(1:end-2)));
% use forward differencing for 1st datapoint to not lose it
Vp(1) = abs((PRES(2)-PRES(1))./(e_time(2)-e_time(1)));
% use backward differencing for last datapoint to not lose it
Vp(end) = abs(((PRES(end)-PRES(end-1))./(e_time(end)-e_time(end-1))));

% convert ascent rate to cm/s
Vp = Vp * 100;

%% Compute thermal mass coefficients based on Vp.
if strcmp(RBRargo_CTcell_model, 'RBRargo|2k_CTcell_2025')
    % These coefficients are for RBRargo|2k CTcell design released in 2025
    ctcoeff = 0.14 * Vp .^ (-1.07);
    alpha   = 0.18 * Vp .^ (-0.54);
    tau     = 37.33 * Vp .^ (-0.50);
elseif strcmp(RBRargo_CTcell_model, 'RBRargo|2k_CTcell_pre2025')
    % These coefficients are for RBRargo|2k CTcell design released in 2017
    ctcoeff = 0.14 * Vp .^ (-1.00);
    alpha   = 0.37 * Vp .^ (-1.03);
    tau     = 16.02 * Vp .^ (-0.26);
elseif strcmp(RBRargo_CTcell_model, 'RBRargo|6k_CTcell')
    % These coefficients are for RBRargo|6k CTcell
    ctcoeff = 0.09 * Vp .^ (-1.17);
    alpha   = 0.27 * Vp .^ (-0.74);
    tau     = 31.54 * Vp .^ (-0.22);
else
    error('RBRargo_CTcell_model must be "RBRargo|2k_CTcell_2025", "RBRargo|2k_CTcell_pre2025" or "RBRargo|6k_CTcell"');
end

%% Compute Tcor (C-T lag)
CTlag = -0.35;
Tcor = interp1(e_time, TEMP, e_time - CTlag);

%% Compute Tlong (long-term thermal mass temperature anomaly)
Tlong = ctcoeff.*(TEMP_CNDC-Tcor);

%% Compute Tshort (long-term thermal mass temperature anomaly)

% Interpolate needed variables onto a 1Hz time series
gooddata = find(~isnan(TEMP) & ~isnan(e_time));
gridded.e_time = min(e_time):1:max(e_time)+1;
gridded.TEMP = interp1(e_time(gooddata),TEMP(gooddata),gridded.e_time);

% Sampling frequency
fs = 1;
% Nyquist frequency
fn = fs/2;

% Lueck and Picklo (1991)'s coefficients (a,b)
a = 4*fn.*alpha.*tau./(1+4*fn.*tau);
b = 1-2.*a./alpha;
gridded.a = interp1(e_time(gooddata),a(gooddata),gridded.e_time);
gridded.b = interp1(e_time(gooddata),b(gooddata),gridded.e_time);

% Apply Lueck and Picklo (1991) recursive filter
gridded.Tshort = zeros(size(gridded.TEMP));

% Set first value to 0°C.
gridded.Tshort(1) = 0;

for tt = 2:length(gridded.TEMP)
    gridded.Tshort(tt) = -gridded.b(tt)*gridded.Tshort(tt-1) + gridded.a(tt)*(gridded.TEMP(tt)-gridded.TEMP(tt-1));
end; clear tt


% Go back and set the first temperature adjustment to the second, since the
% float is probably rising as it takes its first samples
gridded.Tshort(1)=gridded.Tshort(2);

% Re-grid onto original time grid
Tshort = NaN*TEMP;
Tshort(gooddata) = interp1(gridded.e_time,gridded.Tshort,e_time(gooddata));
clear gridded a b fs fn gooddata

%% Compute TEMPcell

% Compute the adjusted temperature that should be used to compute
% salinity corrected for thermal mass errors. It should ONLY be used to
% calculate the adjusted salinity, it does not reflect the true temperature
% of the water. Note the signs, as per the convention established in
% Morison (1994)

TEMPcell = Tcor + Tlong - Tshort;