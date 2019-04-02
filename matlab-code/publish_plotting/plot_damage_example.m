clc
clear

close all
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
% initialize paths.
init_paths();
addpath functions/scanning_v1/
% initialize paths.
init_paths();

Ts = 40e-6;

size_dir = '5microns';

data_root = PATHS.cs_image_data(size_dir, '3-7-2019');

cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-isconnect_out_3-7-2019-01.csv';


chan_map = ChannelMap([1:5]);
% -------------------
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
%
close all
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);

cs_paths = get_cs_paths(data_root, cs_exp_data_name_s{1});

hole_depth = (20);

cs_exp = CsExp(cs_paths, 'feature_height', hole_depth);% 'gg', gg);
cs_exp.uz = detrend(cs_exp.uz);
cs_exp.print_state_times();
%%
figbase = 20;
figure(1); ax1 = gca();
Fig = mkfig(figbase+2, 5, 3); clf
[ha, pos] = tight_subplot(1, 1, [.01], [.15, .02], [.12, .01], false)

cs_exp.meta_exp.z_axis_params.enable_KI_sched

cs_exp.plot_all_cycles(ax1, ha);

ref = cs_exp.meta_exp.z_axis_params.setpoint_scan;
xlm = [3.8566    4.2150];
ylm = [-0.9478    0.4434]
set(ha, 'XLim', xlm, 'YLim', ylm)

% h_ok = ciplot([ylm(1), ylm(1)], [ref, ref], xlm, 'g', ha);
h_no = ciplot([ref, ref], [ylm(2), ylm(2)], xlm, 'r', ha);
% alpha(h_ok, '0.25')
alpha(h_no, '0.25')

% h_ok.DisplayName = 'No Penalty';
h_no.DisplayName = 'Penalize';
xlabel(ha, 'time [s]')
ylabel(ha, 'deflection [V]')
title(ha, '')
%%
leg = legend();
set(leg, 'NumColumns', 2, 'location', 'northwest')

save_fig(Fig, fullfile(PATHS.thesis_fig_final, 'damage_illustration'))
%%
