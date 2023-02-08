function Cnew = compressibility_RBRargo3(COND,PRES,PRES_ADJUSTED,WMO)

%{
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
28 Nov 2021

v1.0 - 28/11/2021
%}

%% checks 

% Check inputs are vectors
[n,m]=size(COND);
if min([n,m])>1
    disp('Vector inputs only, please')
    return
end; clear m n

% Check sizes are right
if size(COND)~=size(PRES) | length(WMO)~=1
    disp('Input size not matching')
    return
end

%% Import look-up table with new coefficients
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 14, "Encoding", "UTF-8");
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["WMO", "DeploymentDate", "Floatmodel", "program", "RBRSN", "X2old", "X3old", "X4old", "X2", "X3", "X4", "Source", "Notes", "VarName14"];
opts.VariableTypes = ["double", "string", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
opts = setvaropts(opts, ["DeploymentDate", "Source", "Notes", "VarName14"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["DeploymentDate", "Floatmodel", "program", "Source", "Notes", "VarName14"], "EmptyFieldRule", "auto");
% Import the data
RBRargo3table = readtable("RBRargo3_compressibility_table.csv", opts);
% Clear temporary variables
clear opts


row = find(RBRargo3table.WMO ==  WMO);
if isempty(row)
    warning('WMO not found in the database')
end

X2old = RBRargo3table.X2old(row);
X3old = RBRargo3table.X3old(row);
X4old = RBRargo3table.X4old(row);

% if new coefficients could not be determined in situ, old coefficients are
% used.
if isnan(RBRargo3table.X2(row)) | isnan(RBRargo3table.X3(row)) | isnan(RBRargo3table.X4(row))
    X2 = X2old;
    X3 = X3old;
    X4 = X4old;
    warning('Sorry, no custom compressibility correction available. Default correction is used!')
else
    X2 = RBRargo3table.X2(row);
    X3 = RBRargo3table.X3(row);
    X4 = RBRargo3table.X4(row);
end


%% Correct conductivity

Cnew = COND .*...
    (1 + X2old*PRES + X3old*PRES.^2 + X4old*PRES.^3)./...
    (1 + X2*PRES_ADJUSTED + X3*PRES_ADJUSTED.^2 + X4*PRES_ADJUSTED.^3);

