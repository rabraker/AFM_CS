clc
clear

close all
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
% initialize paths.
init_paths();
figbase = 20;


size_dir = '5microns';
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
hole_depth = (20);

chan_map = ChannelMap([1:5]);
exp_date1 = '3-21-2019';
exp_date2 = '3-20-2019';
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date1);
cs_file1 = 'cs-traj-512pix-15perc-500nm-5mic-05Hz-250prescan-notconnect_out_3-21-2019-01.csv';
cs_file2 = 'cs-traj-512pix-15perc-500nm-5mic-05Hz-250prescan-notconnect_out_3-20-2019-01.csv';

data_root1 = PATHS.cs_image_data(size_dir, exp_date1);
data_root2 = PATHS.cs_image_data(size_dir, exp_date2);

cs_paths1 = get_cs_paths(data_root1, cs_file1);
cs_paths2 = get_cs_paths(data_root2, cs_file2);

gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);

cs_exps1 = CsExp(cs_paths1, 'feature_height', hole_depth, 'gg', gg);
cs_exps1.print_state_times();

cs_exps2 = CsExp(cs_paths2, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
cs_exps2.print_state_times();


%%
figase=2000;
[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps1.plot_all_cycles(axs{1:4});

%%
% figure(101); ax_phony=gca()
Fig = mkfig(100, 7, 7); clf
ha = tight_subplot(3, 2, [.015, .02], [0.08, .04], [.12, .035], false);

t_start = 8.715;
t_end = 8.755;

indc = {   'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
          'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};

n=24;
CS_idx1 = 274+n;  % cs_exps1.find_cycle_idx(t_start);
CS_idx2 = 275+n;  % cs_exps1.find_cycle_idx(t_end);        
% The Bad Case
traj_names = {'ze', 'y'};
state_names = {'move', 'tdown', 'tsettle', 'scan', 'tup'};

% The Good Case
traj_names = {'uz', 'y', 'x'};
names = {'deflection', 'Y', 'X',};
xlm1 = [9.808, 9.847];
for j=1:length(traj_names)
  h_j = ha( (j-1)*2 + 1 );
  for k=1:length(state_names)
    state_name = state_names{k};
    idx_move = cs_exps1.plot_traj_from_csidx_by_state(CS_idx1, CS_idx2,...
      state_name, traj_names{j}, h_j, 0, 'color', indc{1, k});
  end
  xlim(h_j, xlm1);
  grid(h_j, 'on')
  ylabel(h_j, [names{j} ' [v]'], 'FontSize', 12);
  if j<3
    set(h_j, 'XTickLabel', [])
  else
    xlabel(h_j, 'time [s]')
  end
end

% The Good Case

idx = cs_exps2.idx_state_s.(state_name){CS_idx1};

xlm2 = [11.185, 11.23];
hands = gobjects(length(state_names), 1);
for j=1:length(traj_names)
  h_j = ha( (j-1)*2+2 );
  for k=1:length(state_names)
    state_name = state_names{k};
    hands(k) = cs_exps2.plot_traj_from_csidx_by_state(CS_idx1, CS_idx2,...
      state_name, traj_names{j}, h_j, 0, 'color', indc{1, k});
    xlim(h_j, xlm2)
    hands(k).DisplayName = indc{2, k};
    grid(h_j, 'on')
    if j<3
      set(h_j, 'XTickLabel', [])
    else
      xlabel(h_j, 'time [s]')
    end
  end
end
for k=1:3
   set(ha((k-1)*2+2), 'YTickLabel', '')
end

for k=1:6
   set(ha(k), 'FontSize', 14) 
end

title(ha(1), 'constant-$\rho$')
title(ha(2), 'choose-$\zeta$')

leg = legend(hands);
set(leg, 'location', 'northeast')

%%


save_fig(Fig, fullfile(PATHS.defense_fig(), 'aggresiveD_y_perturbs_z'), true)

%%


