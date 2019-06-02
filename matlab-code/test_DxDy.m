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
exp_date = '5-31-2019'
% ----------------------- Load and Process CS-data -----------------------------
%%
cs_files = {...
%'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_6-1-2019-01.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_6-1-2019-06.csv',...
};

% data_root = '/media/labserver/afm-cs/z-bounce/5-31-2019/';
data_root = 'Z:\afm-cs\z-bounce\6-1-2019';


use_ze = false;
% data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg,...
      'load_full', true, 'reload_raw', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end
%
if 0
[~, axs] = make_cs_traj_figs(1, 2);
cs_exps{1}.plot_all_cycles([], [], axs{:}, [], 512)
%
[~, axs2] = make_cs_traj_figs(102, 2);
cs_exps{2}.plot_all_cycles([], [], axs2{:}, [], 512)
end
X=0;
for k=1:10
idx = cs_exps{1}.idx_state_s.scan{9*k};
% idx = idx-10:idx(end);
t = AFM.Ts*idx;
figure(10); 
xx = cs_exps{1}.x(idx);
plot(t, xx)

figure(11)
plot(t, detrend(xx)*512)

[X_, f] = power_spectrum(detrend(xx), AFM.Ts);
X = X + X_;
end

figure(12)
semilogx(f, 10*log10(X/k))
