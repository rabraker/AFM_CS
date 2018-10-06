clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;
addpath('functions')

% ---------------------------------------------------
% data_root = fullfile(getdataroot(), 'cs-data\20microns\9-23-2017');

data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\cs-data\20microns\parents\';
cs_exp_data_name = 'cs-traj-z-bounce_out_7-11-2018-06.csv';

cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);


dat = csvread(cs_data_path); 

load(cs_meta_path);  % Provides ExpMetaData

fprintf('loading done...\n')

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

close all

x = dat(:,1);
y=dat(:,2);
t = [0:1:length(x)-1]'*Ts;

z_err = dat(:,3);
uz = dat(:,5);
met_ind = dat(:,6);


figbase = 1;
figure(10+figbase)
% plot(x)
title('x')
ax1 = gca;
plotbyindex(ax1, t, x, met_ind, indc);


figure(20+figbase)
% plot(uz)
ax2 = gca
plotbyindex(ax2, t, uz, met_ind, indc);
title('uz')

figure(30+figbase)
% plot(z_err)
ax3 = gca();
plotbyindex(ax3, t, z_err, met_ind, indc);
title('z-err')
hold on
plot([t(1), t(end)], [.05, .05])
plot([t(1), t(end)], -[.05, .05])

linkaxes([ax1, ax2, ax3], 'x')



%%
figure(40+figbase)
ax4=gca();
plotbyindex(ax4, t, y, met_ind, indc);
title('y')
% figure(40)
% plot(min(met_ind, 3))
% title('meta index')
% ax4 = gca();
% linkaxes([ax1, ax2, ax3, ax4], 'x')
linkaxes([ax1, ax2, ax3, ax4], 'x')
