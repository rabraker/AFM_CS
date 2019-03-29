clc
clear
%%
% close all
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

data_root_s = {'/media/labserver/afm-cs/z-bounce/3-28-2019/',...
  '/media/labserver/afm-cs/z-bounce/3-28-2019/'}
%   '/media/labserver/afm-cs/imaging/cs-imaging/5microns/3-28-2019/'}
% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-z-bounce_5000samples_out_3-28-2019-02.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-28-2019-01.csv',...
'cs-traj-512pix-7perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-28-2019-01.csv',...
};

cs_exps = {};
for k=1:2 %length(cs_files)
  cs_paths = get_cs_paths(data_root_s{k}, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  gg = zpk([], [], 1, AFM.Ts);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end


% [~, axs1] = make_cs_traj_figs(figbase, 2);
% cs_exps{1}.plot_all_cycles(axs1{1:2});

[~, axs2] = make_cs_traj_figs(figbase+100, 4);
cs_exps{2}.plot_all_cycles(axs2{1:4});


set(axs1{1}, 'YLim', [-0.25, 0.25])
set(axs1{2}, 'YLim', [-0.4, 0])

set(axs2{1}, 'YLim', [-0.25, 0.25])
set(axs2{2}, 'YLim', [-0.4, 0])
%%
ii_b=1;
N_b = length(cs_exps{ii_b}.idx_state_s.scan);

ZE_b = 0;
for k = 1:N_b
  idxk = cs_exps{ii_b}.idx_state_s.scan{k};
  
  zek = cs_exps{ii_b}.ze(idxk);
  [ZEK,freqs_b] = power_spectrum(zek, AFM.Ts);
  ZE_b = ZEK + ZE_b;
end

ZE_b = ZE_b/N_b;

%%

ii_mu=2;

N_mu = length(cs_exps{ii_mu}.idx_state_s.scan);

ZE_mu = 0;
ZE_mu_skip = 0;
n_skip=0;
for k = 1:N_mu
  idxk = cs_exps{ii_mu}.idx_state_s.scan{k};
  idxk_mv = cs_exps{ii_mu}.idx_state_s.move{k};
  zek = cs_exps{ii_mu}.ze(idxk);
  [ZEK,freqs] = power_spectrum(zek, AFM.Ts);
  ZE_mu = ZEK + ZE_mu;
  
  xk =  cs_exps{ii_mu}.x(idxk_mv);
  if ( max(xk) - min(xk) < 0.2 )
    ZE_mu_skip = ZEK + ZE_mu_skip;
    n_skip = n_skip+1;
  else
    figure(2000)
    plot(xk)
    keyboard
  end
end

ZE_mu = ZE_mu/N_mu;
ZE_mu_skip = ZE_mu_skip/n_skip;

figure(2); clf
semilogx(freqs, 10*log10(abs(ZE_mu)))
hold on
semilogx(freqs, 10*log10(abs(ZE_mu_skip)))

semilogx(freqs_b, 10*log10(abs(ZE_b)))
%%
cs_exp1 = cs_exps{2};
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
