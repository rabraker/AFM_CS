clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;
addpath('functions')
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
data_root = fullfile(getdataroot(), 'cs-data');
data_root = '/media/labserver/acc2018-data/cs-data/5microns/9-22-2017';
cs_exp_data_name = 'cs-traj-512pix-5perc-500nm-5mic-01Hz_out_9-23-2017-04.csv';

% ---------------------------------------------------
% cs_exp_data_name = 'cs-traj10-500_8-22-2017_07.csv';
% cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-06.csv';
% cs_exp_data_name = 'cs-traj10-500_out_8-23-2017-06-meta.csv';



cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);

dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData
%%
% Slice indeces
clc
% k1 = 245210 + 5510;
k1 = 1.9787e+05-200;
k2 =  2.0772e+05;
   
   
% k1 = 250100;
% k2 = 267446;

dat_slice = dat(k1:k2,:);
t = [0:1:length(dat_slice)-1]'*Ts;

x = dat_slice(:,1);
y = dat_slice(:, 2);
z_err = dat_slice(:,3);
uz = dat_slice(:,5);
met_ind = dat_slice(:,6);

% --------------------------- Build Figure ------------------------------
figwidth = 3.4;
figheight = 5.5;
F1 = figure(1); clf;
set(F1, 'Units', 'Inches', 'Position', [-10,0, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');

% indc = {'k', 'r', [0, .75, .75], 'm', [.93 .69 .13], 'b'};
indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};
lft = .157;

ht1 = .24;
ht2 = .37;

ypad = .03;
bt3 = .08;
bt2 = bt3+ht1+ypad;
bt1 = bt2+ht2+ypad;
wd = 1-lft - .025;


ax1 = axes('Position', [lft, bt1, wd, ht1]);
ax2 = axes('Position', [lft, bt2, wd, ht2]);
ax3 = axes('Position', [lft, bt3, wd, ht1]);


F1.CurrentAxes = ax1;

h = plotbyindex(ax1, t, x*volts2microns, met_ind, indc);
h = h(1:5);


ylabel('x [$\mu$m]', 'interpreter', 'latex', 'FontSize', 12)
xlim([t(1), t(end)])

set(ax1, 'XTickLabel', [])
h(4).DisplayName = 'tip engage';
L1 = [h(2), h(3),h(4), h(1)];
% leg1 = legend([h(1), h(2), h(3), h(4), h(5)]);
leg1 = legend(L1);
leg1.Position = [0.6525    0.7795    0.3158    0.1195];
leg1.Interpreter = 'latex';
grid on;


F1.CurrentAxes = ax2;

plotbyindex(ax2, t, uz, met_ind, indc);
ylabel('$u_z$ [v]', 'interpreter', 'latex', 'FontSize', 12)
xlim([t(1), t(end)]);
set(ax2, 'XTickLabel', []);
grid on
ylm = ylim;
ylim([-.2, ylm(2)])


F1.CurrentAxes = ax3;
plotbyindex(ax3, t, z_err, met_ind, indc);
ylabel('deflection [v]', 'interpreter', 'latex', 'FontSize', 12)
xlim([t(1), t(end)])
xlabel('time [s]', 'interpreter', 'latex', 'FontSize', 12)
grid on

linkaxes([ax1, ax2, ax3], 'x')


%%
figname = 'three_CS_cycle.eps'
figpath = fullfile(getfigroot(), figname);

% export_fig(F1, figpath, '-q101');
saveEps(F1, figpath)