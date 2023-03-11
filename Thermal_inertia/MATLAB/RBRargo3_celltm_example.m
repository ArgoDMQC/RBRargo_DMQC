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
argo = readtable("../../test_data/6903078_testdata.csv", opts);

% Clear temporary variables
clear opts

% Simplify variables
lon = argo.LONGITUDEdegree_east(1);
lat = argo.LATITUDEdegree_north(1);
PRES = argo.PRESdecibar;
TEMP = argo.TEMPdegree_Celsius;
PSAL = argo.PSALpsu;
TEMP_CNDC = argo.TEMP_CNDCdegree_Celsius;

%% Infers elptime using a nominal ascent rate of 10 cm/s
elptime = (max(PRES) - PRES)/0.10;
% Sort the data chronologically
[elptime,I] = sort(elptime);
PRES = PRES(I);
TEMP = TEMP(I);
PSAL = PSAL(I);
TEMP_CNDC = TEMP_CNDC(I);

%% Compute additional variables using GSW TEOS-10 
COND = gsw_C_from_SP(PSAL,TEMP,PRES);
SA = gsw_SA_from_SP(PSAL,PRES,lon,lat);
CT = gsw_CT_from_t(SA,TEMP,PRES);
sig = gsw_sigma0(SA,CT);

%%
TEMP_celltm = RBRargo3_celltm(TEMP,PRES,TEMP_CNDC,elptime);

PSAL_ADJUSTED_CTM = gsw_SP_from_C(COND, TEMP_celltm, PRES);

%% Compute additional variables using GSW TEOS-10 
SAcor = gsw_SA_from_SP(PSAL_ADJUSTED_CTM,PRES,lon,lat);
CTcor = gsw_CT_from_t(SAcor,TEMP,PRES);
sigcor = gsw_sigma0(SAcor,CTcor);

%% Plot the difference in conductivity before and after the correction as a function of pressure 
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
plot(TEMP,PRES,'k','linewidth',2);
hold on
plot(TEMP_CNDC,PRES,'--k','linewidth',2);
axis ij
grid on; grid minor
ylim([0 80])
ylabel('P [dbar]')
xlabel('T [°C]')
legend('TEMP','TEMP_CNDC', 'interpreter','none','location','southeast')
set(gca,'fontsize',18, 'xaxislocation','top')

subplot(1,3,2)
plot(PSAL,PRES,'k','linewidth',2);
hold on
plot(PSAL_ADJUSTED_Padj_CTM,PRES,'r','linewidth',2);
axis ij
grid on; grid minor
ylim([0 80])
ylabel('P [dbar]')
xlabel('S [ ]')
legend('PSAL','PSAL_ADJUSTED_CTM','interpreter','none','location','northwest')
set(gca,'fontsize',18, 'xaxislocation','top')

subplot(1,3,3)
plot(sig,PRES,'k','linewidth',2);
hold on
plot(sigcor,PRES,'r','linewidth',2);
axis ij
grid on; grid minor
ylim([0 80])
ylabel('P [dbar]')
xlabel('S [ ]')
legend('\sigma_0','\sigma_0 (ADJUSTED CTM)')
set(gca,'fontsize',18, 'xaxislocation','top')




