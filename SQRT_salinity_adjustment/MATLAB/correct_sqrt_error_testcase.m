clear

data = readmatrix('testdata.csv', 'Range', 2);
Sraw = data(:,5);
Sgood = correct_sqrt_error(Sraw);
%%
figure
subplot(1,2,1)
scatter(data(:,4),data(:,5)-data(:,4),50,data(:,3),'.')
grid on; grid minor
xlim([20,40])
ylim([-.5e-3, 2e-3])
xlabel('Salinity')
ylabel('Salinity error')
fsize(16)
cb = colorbar;
ylabel(cb, 'Temperature [°C]');
title('Before correction')
axis square

subplot(1,2,2)
scatter(data(:,4),Sgood-data(:,4),50,data(:,3),'.')
grid on; grid minor
xlim([20,40])
ylim([-.5e-3, 2e-3])
xlabel('Salinity')
ylabel('Salinity error')
cb = colorbar;
ylabel(cb, 'Temperature [°C]');
title('After correction')
axis square
