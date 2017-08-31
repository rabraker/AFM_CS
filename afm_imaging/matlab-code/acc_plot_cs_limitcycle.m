clc
clear
pix = 256;
width = 5;
microns_per_volt = 50/10;
pix_per_volt = (pix/width)*microns_per_volt;
Ts = 40e-6;

% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\';
data_root = fullfile(getdataroot(), 'cs-data');
% ---------------------------------------------------

cs_exp_data_name = 'cs-traj10-500_out_8-25-2017-06.csv'; % k = 270202, 


cs_exp_meta_name = strrep(cs_exp_data_name, '.csv', '-meta.mat');

cs_data_path = fullfile(data_root, cs_exp_data_name);
cs_meta_path = fullfile(data_root, cs_exp_meta_name);

dat = csvread(cs_data_path); 
load(cs_meta_path);  % Provides ExpMetaData
%%
% Slice indeces
clc
close all
saveall = 1;

k1 = 270202, 
k2 = 290000;

dat_slice = dat(k1:k2,:);
t = [0:1:length(dat_slice)-1]'*Ts;

x = dat_slice(:,1);
y = dat_slice(:, 2);
z_err = dat_slice(:,3);
uz = dat_slice(:,5);
met_ind = dat_slice(:,6);
plot(z_err)
%

dz = ExpMetaData.ramp_rate;
s_dz = sprintf('$\\delta_z = %.3g~[v/Ts]$\n $K_{i,z}=%.3f$', dz, ExpMetaData.Ki);

% --------------------------- Build Figure ------------------------------
figwidth = 3.4;
figheight = 3.;
F1 = figure(1); clf;
set(F1, 'Units', 'Inches', 'Position', [14,6, figwidth, figheight],...
    'PaperUnits', 'Inches', 'PaperSize', [figwidth, figheight])
set(F1, 'Color', 'w');

indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};

lft = .157;
ht1 = .34;
ht2 = .47;

ypad = .03;
bt2 = .14;
bt1 = bt2+ht2+ypad;
% bt1 = bt2+ht2+ypad;
wd = 1-lft - .01;


ax1 = axes('Position', [lft, bt1, wd, ht1]);
ax2 = axes('Position', [lft, bt2, wd, ht2]);
% ax3 = axes('Position', [lft, bt3, wd, ht1]);

linkaxes([ax1,ax2], 'x')
F1.CurrentAxes = ax2;

% h = plotbyindex(ax1, t, x*volts2microns, met_ind, indc);

% ylabel('x [$\mu$m]', 'interpreter', 'latex', 'FontSize', 12)
% xlim([t(1), t(end)])

% set(ax1, 'XTickLabel', [])
% leg1 = legend([h(1), h(2), h(3), h(4), h(5)]);
% grid on;


F1.CurrentAxes = ax1;

plotbyindex(ax1, t, uz, met_ind, indc);

grid on
ylabel('$u_z$ [v]', 'interpreter', 'latex', 'FontSize', 12)
set(ax1, 'XTickLabel', []);
ylm = ylim;
ylim([-.2, ylm(2)])
text(.5,.5, s_dz, 'units', 'normalized', 'interpreter', 'latex')


F1.CurrentAxes = ax2;
h = plotbyindex(ax2, t, z_err, met_ind, indc);
ylabel('deflection [v]', 'interpreter', 'latex', 'FontSize', 12)
xlim([t(1), t(end)])
xlabel('time [s]', 'interpreter', 'latex', 'FontSize', 12)
grid on

ylim([-1, 0.5])

% -------------- legend ---------------
h = h(3:7);
% h(1).DisplayName = 'xy move';
% h(2).DisplayName = 'tip down';
% h(3).DisplayName = 'tip settle';
% h(4).DisplayName = '$\mu$-path scan';
% h(5).DisplayName = 'tip up';
L1 = [h(1), h(2), h(3)];
L2 = [h(4), h(5)];

leg1 = legend([L1, L2]);%%
leg1.Position = [0.6197 0.1629 0.3690 0.2708];
leg1.Interpreter = 'latex';



linkaxes([ax1, ax2], 'x')


%%
if saveall
    figname = 'CS_limitcycle.pdf'
    figpath = fullfile(getfigroot(), figname);

    export_fig(F1, figpath, '-q101');

end


