%{
This routine is provided as an example on how to correct conductivity data
collected from an RBRargo3 Argo float for compressibility errors.

It relies on the CSV file listing the appropriate coefficients:
https://github.com/ArgoDMQC/RBRargo_DMQC/blob/main/RBRargo3_compressibility_table.csv

As well as on an example profile from:
https://fleetmonitoring.euro-argo.eu/float/6903078
%}

clear
%% Load the CSV file with example data

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["PLATFORM_CODE", "DATEYYYYMMDDTHHMISSZ", "DATE_QC", "LATITUDEdegree_north", "LONGITUDEdegree_east", "POSITION_QC", "PRESdecibar", "PRES_QC", "PSALpsu", "PSAL_QC", "TEMPdegree_Celsius", "TEMP_QC", "PRES_ADJUSTEDdecibar", "PRES_ADJUSTED_QC", "TEMP_ADJUSTEDdegree_Celsius", "TEMP_ADJUSTED_QC", "PSAL_ADJUSTEDpsu", "PSAL_ADJUSTED_QC", "TEMP_CNDCdegree_Celsius", "TEMP_CNDC_QC"];
opts.VariableTypes = ["double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "DATEYYYYMMDDTHHMISSZ", "EmptyFieldRule", "auto");

% Import the data
argo = readtable("6903078_testdata.csv", opts);

% Clear temporary variables
clear opts

%% Compute conductivity from saliniy using GSW TEOS-10 
argo.COND = gsw_C_from_SP(argo.PSALpsu,argo.TEMPdegree_Celsius,argo.PRESdecibar);

%% Applies the new compressibility correction for that particular float
Cnew = compressibility_RBRargo3(argo.COND,argo.PRESdecibar,argo.PRES_ADJUSTEDdecibar,6903078);

%% Plot the difference in conductivity before and after the correction as a function of pressure 

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
plot(argo.COND-Cnew,argo.PRESdecibar,'k','linewidth',2);
axis ij
grid on; grid minor
ylim([0 2000])
ylabel('P [dbar]')
xlabel('∆C [mS/cm]')
title({'Change in conductivity after applying','custom compressibility correction'})
set(gca,'fontsize',18)

subplot(1,3,2)
plot(argo.PSALpsu-gsw_SP_from_C(Cnew,argo.TEMPdegree_Celsius,argo.PRESdecibar),argo.PRESdecibar,'k','linewidth',2);
axis ij
grid on; grid minor
ylim([0 2000])
ylabel('P [dbar]')
xlabel('∆S [ ]')
title({'Change in salinity after applying','custom compressibility correction'})
set(gca,'fontsize',18)

subplot(1,3,3)
plot(argo.PSALpsu,argo.PRESdecibar,'-k','linewidth',2);
hold on
plot(gsw_SP_from_C(Cnew,argo.TEMPdegree_Celsius,argo.PRESdecibar),argo.PRESdecibar,'-r','linewidth',2);
axis ij
grid on; grid minor
ylim([0 2000])
ylabel('P [dbar]')
xlabel('PSAL [ ]')
%title({'Change in salinity after applying','custom compressibility correction'})
set(gca,'fontsize',18)
legend('PSAL','PSAL_ ADJUSTED','location','southeast')





