%{
This routine is provided as an example on how to infer internal temperature
using the temperature measured from an RBRargo3 Argo.

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

COND = gsw_C_from_SP(argo.PSAL_ADJUSTEDpsu,argo.TEMP_ADJUSTEDdegree_Celsius,argo.PRES_ADJUSTEDdecibar);


%% Infers elptime using a nominal ascent rate of 10 cm/s

elptime = (max(argo.PRES_ADJUSTEDdecibar) - argo.PRES_ADJUSTEDdecibar)/0.10;

% Sort the data chronologically
[elptime,I] = sort(elptime);
TEMP_ADJUSTED = argo.TEMP_ADJUSTEDdegree_Celsius(I);

TEMP_CNDC = RBRargo3_TEMP_CNDC_from_TEMP_ADJUSTED(TEMP_ADJUSTED,elptime);
TEMP_CNDC = flipud(TEMP_CNDC);
%%

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
plot(argo.TEMP_ADJUSTEDdegree_Celsius,argo.PRESdecibar,'k','linewidth',2);
hold on
plot(argo.TEMP_CNDCdegree_Celsius,argo.PRESdecibar,'--k','linewidth',2);
plot(TEMP_CNDC,argo.PRESdecibar,'--r','linewidth',2);
axis ij
grid on; grid minor
ylim([0 2000])
ylabel('P [dbar]')
xlabel('T [°C]')
legend('TEMP','TEMP_CNDC (measured)','TEMP_CNDC (modeled)', 'interpreter','none','location','southeast')
%title({'Change in conductivity after applying','custom compressibility correction'})
set(gca,'fontsize',18, 'xaxislocation','top')

subplot(1,3,2)
plot(argo.TEMP_CNDCdegree_Celsius-TEMP_CNDC,argo.PRESdecibar,'k','linewidth',2);
axis ij
grid on; grid minor
ylim([0 2000])
ylabel('P [dbar]')
xlabel('∆T [°C]')
XLIM = get(gca,'xlim');
XLIMmax = max(abs(XLIM));
xlim([-XLIMmax XLIMmax]);
%title({'Change in salinity after applying','custom compressibility correction'})
set(gca,'fontsize',18, 'xaxislocation','top')

ctcoeff = 1.4e-2;
subplot(1,3,3)
plot((gsw_SP_from_C(COND,argo.TEMP_ADJUSTEDdegree_Celsius+ctcoeff*(argo.TEMP_CNDCdegree_Celsius-argo.TEMP_ADJUSTEDdegree_Celsius),argo.PRES_ADJUSTEDdecibar)-...
    gsw_SP_from_C(COND,argo.TEMP_ADJUSTEDdegree_Celsius+ctcoeff*(TEMP_CNDC-argo.TEMP_ADJUSTEDdegree_Celsius),argo.PRES_ADJUSTEDdecibar)),...
    argo.PRES_ADJUSTEDdecibar,'k','linewidth',2)
grid on
axis ij
axis tight
ylim([0 2000])
ylabel('PRES [dbar]')
xlabel('∆S [ ]')
title('Salinity error using TEMP CNDC (modeled)')
set(gca,'fontsize',18, 'xaxislocation','top')
XLIM = get(gca,'xlim');
XLIMmax = max(abs(XLIM));
xlim([-XLIMmax XLIMmax]);
% plot(argo.PSALpsu,argo.PRESdecibar,'-k','linewidth',2);
% hold on
% plot(gsw_SP_from_C(Cnew,argo.TEMPdegree_Celsius,argo.PRESdecibar),argo.PRESdecibar,'-r','linewidth',2);
% axis ij
% grid on; grid minor
% ylim([0 2000])
% ylabel('P [dbar]')
% xlabel('PSAL [ ]')
% %title({'Change in salinity after applying','custom compressibility correction'})
% set(gca,'fontsize',18)
% legend('PSAL','PSAL_ADJUSTED_Padj','location','southeast','Interpreter','none')
