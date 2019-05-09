clc
clear
%%
% close all
clear PATHS
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
addpath ~/matlab/afm-cs/matlab-code/functions/
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
% exp_date = '3-20-2019'
% ----------------------- Load and Process CS-data -----------------------------

data_root_s = {'/media/labserver/afm-cs/z-bounce/4-12-2019/',...
  '/media/labserver/afm-cs/z-bounce/4-12-2019/'};
%   '/media/labserver/afm-cs/imaging/cs-imaging/5microns/3-28-2019/'}
% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-z-bounce_1500samples_Npath160_prescan250_out_4-12-2019-01.csv',...
'cs-traj-512pix-19perc-500nm-5mic-01Hz-250prescan-notconnect_out_4-12-2019-01.csv',...
};

cs_exps = {};
for k=1:2 %length(cs_files)
  cs_paths = get_cs_paths(data_root_s{k}, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
%   gg = zpk([], [], 1, AFM.Ts);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

%%
[~, axs1] = make_cs_traj_figs(figbase, 2);
cs_exps{1}.plot_all_cycles(axs1{1:2});

[~, axs2] = make_cs_traj_figs(figbase+100, 4);
cs_exps{2}.plot_all_cycles(axs2{1:4}, [], 512);

linkaxes([axs1{:}, axs2{:}], 'x')
%%
% set(axs1{1}, 'YLim', [-0.25, 0.25])
% set(axs1{2}, 'YLim', [-0.4, 0])
% 
% set(axs2{1}, 'YLim', [-0.25, 0.25])
% set(axs2{2}, 'YLim', [-0.4, 0])

for k=1:length(cs_exps)
   cs_exps{k}.uz  = cs_exps{k}.uz*AFM.volts2nm_z;
end
ii_b=1;
N_b = length(cs_exps{ii_b}.idx_state_s.scan);

UZ_b = 0;
for k = 1:N_b
  idxk = cs_exps{ii_b}.idx_state_s.scan{k};
  
  uzk = cs_exps{ii_b}.uz(idxk);
  [UZK,freqs_b] = power_spectrum_local(uzk); %, AFM.Ts);
  UZ_b = UZK + UZ_b;
end

UZ_b = UZ_b/N_b;



ii_mu=2;

N_mu = length(cs_exps{ii_mu}.idx_state_s.scan);

UZ_mu = 0;
UZ_mu_skip = 0;
UZ_mu_jump = 0;

XK_skip=0;
YK_skip=0;

n_skip=0;
n_jump=0;
skips = [9:8:160*5];
% for k = 1:N_mu
%   idxk = cs_exps{ii_mu}.idx_state_s.scan{k};
%   idxk_mv = cs_exps{ii_mu}.idx_state_s.move{k};
%   uzk = cs_exps{ii_mu}.uz(idxk);
%   [UZK,freqs] = power_spectrum_local(uzk); %, AFM.Ts);
%   UZ_mu = UZK + UZ_mu;
%   
%   xk =  cs_exps{ii_mu}.x(idxk_mv);
%   
%   if isempty(intersect(k, skips)) && isempty(intersect(k, skips+1))
%     UZ_mu_skip = UZK + UZ_mu_skip;
%     n_skip = n_skip+1;
%   else
%     figure(2000)
%     plot(xk)
%     k
% %     keyboard
%   end
% end

for k = 1:N_mu
  idxk = cs_exps{ii_mu}.idx_state_s.scan{k};
  idxk_mv = cs_exps{ii_mu}.idx_state_s.move{k};
  uzk = cs_exps{ii_mu}.uz(idxk);
  xk = cs_exps{ii_mu}.x(idxk);
  yk = cs_exps{ii_mu}.y(idxk);
  [UZK,freqs] = power_spectrum_local(uzk); %, AFM.Ts);
  [XK,freqs] = power_spectrum_local(xk); %, AFM.Ts);
  [YK,freqs] = power_spectrum_local(yk); %, AFM.Ts);
  
  UZ_mu = UZK + UZ_mu;
  
  xk =  cs_exps{ii_mu}.x(idxk_mv);
  
  if ~isempty(intersect(k, skips-1))
    UZ_mu_skip = UZK + UZ_mu_skip;
    XK_skip = XK + XK_skip;
    YK_skip = YK + YK_skip;
    n_skip = n_skip + 1;
  end
  
  if ~isempty(intersect(k, skips))
        UZ_mu_jump = UZK + UZ_mu_jump;
        n_jump = n_jump + 1;
  end
end


UZ_mu = UZ_mu/N_mu;
UZ_mu_skip = UZ_mu_skip/n_skip;
UZ_mu_jump = UZ_mu_jump/n_jump;


width = 3.5;
height = 3;
F3 = mkfig(23, width, height); clf
gap = 0.1;
margh = [0.13, .025];
margw = [0.15, .025];

ax3 = tight_subplot(1, 1, gap, margh, margw, false);

F2 = mkfig(20, width, height); clf
ax2 = tight_subplot(1, 1, gap, margh, margw, false);

h0 = semilogx(ax2, freqs, 10*log10(abs(UZ_mu_skip)));
h0.DisplayName = 'no-jump'
hold(ax2, 'on')
h1 = semilogx(ax2, freqs, 10*log10(abs(UZ_mu_jump)));
h1.DisplayName = '5-micron jump'


% h1 = semilogx(freqs, 10*log10(abs(UZ_mu)));
% h1.DisplayName = '$\mu$-path (whole)'
% hold on
h2 = semilogx(ax3, freqs, 10*log10(abs(UZ_mu_skip)));
hold(ax3, 'on')
% semilogx(ax3, freqs, 10*log10(XK_skip))
% semilogx(ax3, freqs, 10*log10(YK_skip))

h2.DisplayName = '$\mu$-path (pre-scan+scan)'


h3 = semilogx(ax3, freqs_b, 10*log10(abs(UZ_b)));
h3.DisplayName = 'z-bounce (no scan)'
leg2 = legend([h2, h3], 'Position', [0.1481 0.1288 0.7262 0.1712], 'FontSize', 14);
grid(ax3, 'on')
xlabel(ax3, 'Hz')
ylabel(ax3, 'PSD')

grid(ax2, 'on')
xlabel(ax2, 'Hz')
ylabel(ax2, 'PSD')
leg1 = legend([h0, h1], 'Position', [0.1492 0.1286 0.5084 0.1712], 'FontSize', 14);
%%
save_fig(F3, fullfile(PATHS.defense_fig(), 'prescan_zbounce_PSD'), true)
save_fig(F3, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/prescan_zbounce_PSD'), false)

save_fig(F2, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/prescan_5mic_vs_still_PSD'), false)
%%
width = 3.5;
height = 3;
gap = 0.1;
margh = 0.1;
margw = 0.1;
margh = [0.15, .1];
margw = [0.15, .04];

F1 = mkfig(21, width, height); clf
ax1 = tight_subplot(1, 1, gap, margh, margw, false);

F2 = mkfig(22, width, height); clf
ax2 = tight_subplot(1, 1, gap, margh, margw, false);
%
to = 1;
tstart = 4.64; %3.46 + to;
tend = 4.68; %3.5 + to;

zbounce = true;
plot_uz_time_range(cs_exps{1}, ax1, tstart, tend, zbounce, 'uz');
zbounce = false;
h = plot_uz_time_range(cs_exps{2}, ax2, tstart, tend, zbounce, 'uz');

grid(ax1, 'on')
grid(ax2, 'on')
ylim(ax1, [-5, 30])
ylim(ax2, [-5, 30])
title(ax1, 'z-bounce (no scan)', 'FontSize', 14)
title(ax2, '$\mu$-path (pre-scan + scan)', 'FontSize', 14)

% delt = 
xlim(ax1, [tstart, tend+.04])
xlim(ax2, [tstart-.04, tend])
s1 = text(ax2, 4.6224, 7.0357, 'pre-scan', 'FontSize', 14)
an1 = annotation('arrow', [0.4405 0.3810], [0.3750 0.2674])

ylabel(ax1, '$u_Z$ [nm]', 'FontSize', 14)
ylabel(ax2, '$u_Z$ [nm]', 'FontSize', 14)
xlabel(ax1, 'time [s]', 'FontSize', 14)
xlabel(ax2, 'time [s]', 'FontSize', 14)

leg = legend(h, 'Position', [0.4921 0.5415 0.3977 0.3542], 'FontSize', 12);
%%
save_fig(F1, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/zbounce_uz_example'), false)
save_fig(F2, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/prescan_uz_example'), false)

save_fig(F2, fullfile(PATHS.defense_fig(), 'prescan_uz_example'), true)
save_fig(F1, fullfile(PATHS.defense_fig(), 'zbounce_uz_example'), true)

%%
title(ax2, 'starting height = 105 nm')
s1.Visible = 'off';
an1.Visible = 'off';

save_fig(F2, fullfile(PATHS.defense_fig(), 'Zdir_pzcancel_works'), true)

%%

F5 = mkfig(25, width, height); clf
axs2 = tight_subplot(2, 1, gap, margh, margw, false);

% title(axs1(1), '$xy$ fixed', 'FontSize', 14)

zbounce = false;
plot_uz_time_range(cs_exps{2}, axs2(1), tstart, tend, zbounce, 'uz');
zbounce = false;
h = plot_uz_time_range(cs_exps{2}, axs2(2), tstart, tend, zbounce, 'x');


grid(axs2(1), 'on')
grid(axs2(2), 'on')

ylim(axs2(1), [-5, 20])


xlim(axs2(1), [4.61, 4.78])
xlim(axs2(2), [4.61, 4.78])

ylim(axs2(2), [-0.5, 5])

ylabel(axs2(1), '$u_Z$~[nm]', 'FontSize', 14)
ylabel(axs2(2), '$x$~[$\mu$m]', 'FontSize', 14)
xlabel(axs2(2), 'time [s]', 'FontSize', 14)

legend(h, 'location', 'northeast');
%%

save_fig(F5, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/prescan_uz_mv_example'), false)
save_fig(F5, fullfile(PATHS.defense_fig(), 'prescan_uz_mv_example'), true)
%%
function hands = plot_uz_time_range(cse, ax, tstart, tend, zbounce, signal)
indc = {'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
        'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};
    if signal == 'x'
        scl = 5;
    else
        scl = 1;
    end
    hold(ax, 'on')
    state_s = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
    hands = gobjects(length(state_s), 1);
    for k=1:length(state_s)
        if k==4 && zbounce
            clr = indc{1,3};
        else
            clr = indc{1,k};
        end
        state = state_s{k};
        
        idxs = cse.get_idx_by_state_in_time_range(state, tstart, tend);
        
        for j=1:length(idxs)
            t = cse.t(idxs{j});
            uzk = cse.(signal)(idxs{j});
            h = plot(ax, t, uzk*scl, 'Color', clr);
        end
        h.DisplayName = indc{2, k};
        hands(k) = h;
    end
end


function [Pxx, freqs] = power_spectrum_local(X)

  N = length(X);

  [Pxx, freqs] = periodogram(detrend(X), [], N, 1/AFM.Ts);
end



