clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;

% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
% data_root = fullfile(getdataroot(), 'cs-data');
data_root = fullfile(getdataroot(), 'cs-data');
% ---------------------------------------------------
% cs_exp_data_name = 'cs-traj10-500_8-22-2017_07.csv';
% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-06.csv';

% Shows a few large oscillations, not quite limit cycle
% cs_exp_data_name = 'cs-traj10-500_out_8-24-2017-04.csv';

% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-06.csv'; % k = 270202, 

% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-09.csv'; % k =
cs_exp_data_name = 'cs-traj-8perc-500nm-5mic-1Hz_out_9-3-2017-02.csv';
% cs_exp_data_name = 'data-in_single_out_9-4-2017-01.csv';
cs_exp_data_name = 'cs-traj-8perc-500nm-5mic-1Hz_out_9-4-2017-04.csv'
cs_exp_data_name = 'cs-traj-10perc-500nm-5mic-1Hz_out_9-4-2017-02.csv';
cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);

dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

%
%

x = dat(:,1);
t = [0:1:length(x)-1]'*Ts;

z_err = dat(:,3);
uz = dat(:,5);
met_ind = dat(:,6);



figure(10)
% plot(x)
title('x')
ax1 = gca;
plotbyindex(ax1, t, x, met_ind, indc);


figure(20)
% plot(uz)
ax2 = gca
plotbyindex(ax2, t, uz, met_ind, indc);
title('uz')

figure(30)
% plot(z_err)
ax3 = gca();
plotbyindex(ax3, t, z_err, met_ind, indc);
title('z-err')

% figure(40)
% plot(min(met_ind, 3))
% title('meta index')
% ax4 = gca();
% linkaxes([ax1, ax2, ax3, ax4], 'x')
linkaxes([ax1, ax2, ax3], 'x')

%%
Ki = 0.005;
Dz = zpk([0], [1], Ki, 40e-6)

frfdat = csvread('C:\Users\arnold\Desktop\rnoise_z_.csv', 1);
frfdat = frfdat(2:end,:);

magDB = frfdat(:,2);
mag = 10.^(magDB/20);

phase_deg = frfdat(:,4);
phase_rad = phase_deg*pi/180;

freq_hz = frfdat(:,1);
freq_rad = freq_hz/2/pi;
g_frf = mag.*exp(j*phase_rad);
% -----------
dz_frf = squeeze(freqresp(Dz, freq_rad));

h_frf = dz_frf.*g_frf./(1 + dz_frf.*g_frf);

figure(50); clf
% subplot(2,1,1)
semilogx(freq_hz, 20*log10(abs(g_frf)))
hold on
semilogx(freq_hz, 20*log10(abs(h_frf)))
semilogx(freq_hz, 20*log10(abs(dz_frf)));
% subplot(2,1,2)
% semilogx(freq_hz, phase_deg)












