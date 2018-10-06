clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;
addpath('functions')
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
% data_root = fullfile(getdataroot(), 'cs-data');
data_root = fullfile(getdataroot(), 'cs-data\20microns\9-23-2017');
data_root = '/media/labserver'
% ---------------------------------------------------
% cs_exp_data_name = 'cs-traj10-500_8-22-2017_07.csv';
% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-06.csv';
% Shows a few large oscillations, not quite limit cycle
% cs_exp_data_name = 'cs-traj10-500_out_8-24-2017-04.csv';
% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-06.csv'; % k = 270202, 
% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-09.csv'; % k =
% cs_exp_data_name = 'cs-traj-8perc-500nm-5mic-1Hz_out_9-3-2017-02.csv';
% cs_exp_data_name = 'data-in_single_out_9-4-2017-01.csv';
% cs_exp_data_name = 'cs-traj10-500_8-22-2017_05.csv';

% data_root = '/media/labserver/acc2018-data/cs-data/5microns/9-22-2017';
% cs_exp_data_name = 'cs-traj-512pix-5perc-500nm-5mic-01Hz_out_9-23-2017-06.csv';
% cs_exp_data_name = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_9-23-2017-04.csv'; % the 10% I used
cs_exp_data_name = 'cs-traj-z-bounce_out_10-5-2018-46.csv';

cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);


dat = csvread(cs_data_path); 

load(cs_meta_path);  % Provides ExpMetaData

fprintf('loading done...\n')

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

%
%%
% close all

x = dat(:,1);
y=dat(:,2);
t = [0:1:length(x)-1]'*Ts;

z_err = dat(:,3);
uz = dat(:,5);
met_ind = dat(:,6);


figbase = 1;
figure(10+figbase), clf
% plot(x)
title('x')
ax1 = gca;
plotbyindex(ax1, t, x, met_ind, indc);


figure(20+figbase), clf
% plot(uz)
ax2 = gca
plotbyindex(ax2, t, uz, met_ind, indc);
title('uz')

figure(30+figbase), clf
% plot(z_err)
ax3 = gca();
hp = plotbyindex(ax3, t, z_err, met_ind, indc);
title('z-err')
hold on
plot([t(1), t(end)], [.05, .05])
plot([t(1), t(end)], -[.05, .05])
legend(hp(1:5))
linkaxes([ax1, ax2, ax3], 'x')




figure(40+figbase), clf
ax4=gca();
plotbyindex(ax4, t, y, met_ind, indc);
title('y')
% figure(40)
% plot(min(met_ind, 3))
% title('meta index')
% ax4 = gca();
% linkaxes([ax1, ax2, ax3, ax4], 'x')
linkaxes([ax1, ax2, ax3, ax4], 'x')

%%
clc
Ts = 40e-6;
tdown_idx = find(met_ind == -5);

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

met_ind0 = met_ind;
met_ind0(met_ind0 >0) = -6;

figure(31)
yyaxis right
plot(t, met_ind0, 'g')

%%

figure(31)

plot(t(tdown_idx), z_err(tdown_idx))

dt = diff(t(tdown_idx));
jump_idx = find(dt > 2*Ts)

typ = -5;
met_ind_tmp = met_ind;
zerr_s = {};
t_s = {};
%          -1       -2       -3       -6       -5
names = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
idx_state_s.scan = {}
idx_state_s.move = {}
idx_state_s.tdown = {}
idx_state_s.tsettle = {}
idx_state_s.tup = {}

idx
for k=1:5
  
  
  
  
end

hold on
for k=1:5
  plot(t_s{k}, zerr_s{k});
  
end



