% This script was supposed to plot the stuff about dynamic detrending. It
% doesnt show anything useful (yet?).

clc
clear

% close all
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
% initialize paths.
init_paths();
addpath functions/scanning_v1/
% initialize paths.
init_paths();

fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);

size_dir = '5microns';
data_root = PATHS.cs_image_data(size_dir, '3-12-2019');

cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-isconnect_out_3-12-2019-01.csv';
chan_map = ChannelMap([1:5]);
% -------------------

cs_paths = get_cs_paths(data_root, cs_exp_data_name_s{1});

hole_depth = (20);

cs_exp1 = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', @log_creep_detrend); %, 'gg', gg);
cs_exp2 = CsExp(cs_paths, 'feature_height', hole_depth);
% cs_exp.uz = detrend(cs_exp.uz);
cs_exp1.print_state_times();
%%
figbase = 20;
figure(1+figbase);clf; 
ax1 = gca();
Fig = mkfig(figbase+2, 7, 4); clf
[ha, pos] = tight_subplot(1, 1, [.01], [.15, .02], [.12, .01], false);


Fig3 = mkfig(figbase+3, 7, 4); clf
[ha3, pos] = tight_subplot(1, 1, [.01], [.15, .02], [.12, .01], false);
figure(4+figbase);clf; 
ax4 = gca();



cs_exp1.plot_all_cycles(ha, ax1); %, ax3, ax4);
cs_exp2.plot_all_cycles(ha3, ax1); %, ax3, ax4);

ref = cs_exp1.meta_exp.z_axis_params.setpoint_scan;
xlm = [19.7876,   21.0337];
ylm = [-0.1418,    0.2663];

% set(ha, 'XLim', xlm, 'YLim', ylm)


xlabel(ha, 'time [s]')
ylabel(ha, '$u_{z_I} [V]$')
title(ha, '')
%%
leg = legend();
set(leg, 'NumColumns', 2, 'location', 'northwest')

save_fig(Fig, 'notes/figures/dinv_dynamic_detrend')
%%
t_ = (0:length(cs_exp1.uz)-1)'*AFM.Ts;
t0 = 0.6940;
t1 = 0.7670;



idx1 = find(t_ > t0, 1, 'first');
idx2 = find(t_ < t1, 1, 'last');

uz = cs_exp1.uz(idx1:idx2);
uz = -(uz - uz(1));
t = t_(idx1:idx2);
ze = cs_exp1.ze(idx1:idx2);
ze = (ze - ze(1));

figure(10)
subplot(2,1,1)
plot(t, uz)
grid on
subplot(2,1,2)
plot(t, ze)



z0=[0.9939, 0.9993];
p0 = [0.9940, 0.9994];
k0= [1]

theta0 = [z0, p0, k0];
gv = zpk([], [], 1, AFM.Ts);
fun = @(theta) fit_gdrift(theta, gv, ze, uz, t, 2);

theta=lsqnonlin(fun, theta0);
num=theta(1:2);
z = theta(1:2);
p = theta(3:4);
k = theta(5);
gg_fit = zpk(z, p, k, AFM.Ts);


y = lsim(gg_fit, uz, t);

subplot(2,1,2)
hold on
plot(t,y )











