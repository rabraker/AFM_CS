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

cs_exp = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', @log_creep_detrend); %, 'gg', gg);
cs_exp.uz = detrend(cs_exp.uz);
cs_exp.print_state_times();
%%
figbase = 20;
figure(1+figbase);clf; 
ax1 = gca();


figure(3+figbase);clf; 
ax3 = gca();
figure(4+figbase);clf; 
ax4 = gca();

Fig = mkfig(figbase+2, 7, 4); clf
[ha, pos] = tight_subplot(1, 1, [.01], [.15, .02], [.12, .01], false);

cs_exp.meta_exp.z_axis_params.enable_KI_sched

cs_exp.plot_all_cycles(ha, ax1); %, ax3, ax4);

ref = cs_exp.meta_exp.z_axis_params.setpoint_scan;
xlm = [19.7876,   21.0337];
ylm = [-0.1418,    0.2663];

set(ha, 'XLim', xlm, 'YLim', ylm)


xlabel(ha, 'time [s]')
ylabel(ha, '$u_{z_I} [V]$')
title(ha, '')
%%
leg = legend();
set(leg, 'NumColumns', 2, 'location', 'northwest')

save_fig(Fig, 'notes/figures/dinv_CS_cycle_decay')
%%
