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

data_root_s = {'/media/labserver/afm-cs/z-bounce/3-29-2019/',...
  '/media/labserver/afm-cs/z-bounce/3-29-2019/'};
%   '/media/labserver/afm-cs/imaging/cs-imaging/5microns/3-28-2019/'}
% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-z-bounce_1500samples_Npath160_prescan250_out_3-29-2019-xy-off-01.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-29-2019-01.csv',...
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


[~, axs1] = make_cs_traj_figs(figbase, 2);
cs_exps{1}.plot_all_cycles(axs1{1:2});

[~, axs2] = make_cs_traj_figs(figbase+100, 4);
cs_exps{2}.plot_all_cycles(axs2{1:4});
% 
% 
% set(axs1{1}, 'YLim', [-0.25, 0.25])
% set(axs1{2}, 'YLim', [-0.4, 0])
% 
% set(axs2{1}, 'YLim', [-0.25, 0.25])
% set(axs2{2}, 'YLim', [-0.4, 0])
%%
ii_b=1;
N_b = length(cs_exps{ii_b}.idx_state_s.scan);

ZE_b = 0;
for k = 1:N_b
  idxk = cs_exps{ii_b}.idx_state_s.scan{k};
  
  zek = cs_exps{ii_b}.uz(idxk);
  [ZEK,freqs_b] = power_spectrum_local(zek); %, AFM.Ts);
  ZE_b = ZEK + ZE_b;
end

ZE_b = ZE_b/N_b;



ii_mu=2;

N_mu = length(cs_exps{ii_mu}.idx_state_s.scan);

ZE_mu = 0;
ZE_mu_skip = 0;
n_skip=0;
skips = [9:8:160];
for k = 1:N_mu
  idxk = cs_exps{ii_mu}.idx_state_s.scan{k};
  idxk_mv = cs_exps{ii_mu}.idx_state_s.move{k};
  zek = cs_exps{ii_mu}.uz(idxk);
  [ZEK,freqs] = power_spectrum_local(zek); %, AFM.Ts);
  ZE_mu = ZEK + ZE_mu;
  
  xk =  cs_exps{ii_mu}.x(idxk_mv);
  
  if isempty(intersect(k, skips)) && isempty(intersect(k, skips+1))
    ZE_mu_skip = ZEK + ZE_mu_skip;
    n_skip = n_skip+1;
  else
    figure(2000)
    plot(xk)
    k
%     keyboard
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



function [Pxx, freqs] = power_spectrum_local(X)

  N = length(X);

  [Pxx, freqs] = periodogram(detrend(X), [], N, 1/AFM.Ts);
end



