clc
clear
%%
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
exp_date = '3-20-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-04Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-05Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-02Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-04Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-05Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
};


data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

%%
[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps{1}.plot_all_cycles(axs{1:4});



%%
cs_exp1 = cs_exps{1};
indc = {   'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
          'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};
        
Fig = mkfig(100, 7, 7); clf
ha = tight_subplot(4, 1, [.05, .05], [0.075, .03], [.08, .01], false);



CS_idx1 = cs_exps{1}.find_cycle_idx(16.31);
CS_idx2 = CS_idx1 + 3;

%
% The Bad Case
traj_names = {'ze', 'y'};
state_names = {'move', 'tdown', 'tsettle', 'scan', 'tup'};

% The Good Case
traj_names = {'y', 'x', 'uz', 'ze'};
names = {'$Y$', '$X$', '$u_Z$', 'deflection'};
xlm1 = [16.25, 16.5];
for j=1:length(traj_names)
  h_j = ha(j);
  for k=1:length(state_names)
    state_name = state_names{k};
    idx_move = cs_exp1.plot_traj_from_csidx_by_state(CS_idx1, CS_idx2,...
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

%%
save_fig(Fig, 'notes/figures/cs_cycle')

save_fig(Fig, fullfile(PATHS.thesis_root, 'plots-afm-cs-final/figures/cs_cycle'))
