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

data_root = '/media/labserver/afm-cs/z-bounce/5-28-2019/';

% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-28-2019-01.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-28-2019-02.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-28-2019-03.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-28-2019-04.csv',...
};

cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
%   gg = zpk([], [], 1, AFM.Ts);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

name_s = {'CR w/o notch', 'CR-w/notch', 'CZ w/o notch', 'CZ w/ notch'};
for k=1:length(cs_exps)
   cs_exps{k}.uz  = cs_exps{k}.uz*AFM.volts2nm_z;
end

%%
[~, axs1] = make_cs_traj_figs(figbase+100, 3);
cs_exps{1}.plot_all_cycles(axs1{1:3}, [], [], 512);

[~, axs2] = make_cs_traj_figs(figbase+200, 3);
cs_exps{2}.plot_all_cycles(axs2{1:3}, [], [], 512);

[~, axs3] = make_cs_traj_figs(figbase+300, 3);
cs_exps{3}.plot_all_cycles(axs3{1:3}, [], [], 512);

[~, axs4] = make_cs_traj_figs(figbase+400, 3);
cs_exps{4}.plot_all_cycles(axs4{1:3}, [], [], 512);

linkaxes([axs1{:}, axs2{:}, axs3{:}, axs4{:}], 'x')
%%
% set(axs1{1}, 'YLim', [-0.25, 0.25])
% set(axs1{2}, 'YLim', [-0.4, 0])
% 
% set(axs2{1}, 'YLim', [-0.25, 0.25])
% set(axs2{2}, 'YLim', [-0.4, 0])


%%
UZ_mu = 0;
UZ_mu_skip = 0;
UZ_mu_jump = 0;

XK_skip=0;
YK_skip=0;

n_skip=0;
n_jump=0;
skips = [9:8:160*5];

UZ_s = {}; 
UZ_pre_jump_s = {}; 
UZ_jump_s = {};
for k=1:length(cs_exps)
    [UZ_, UZ_pre_jump_, UZ_jump_, freqs] = psd_around_jumps(cs_exps{k}, skips);
    UZ_s{k} = UZ_;
    UZ_pre_jump_s{k} = UZ_pre_jump_;
    UZ_jump_s{k} = UZ_jump_;

end
%

figure(1); clf;
subplot(1,2,1)
ax1 = gca();
subplot(1,2,2);
ax2 = gca();

hands_1 = gobjects(length(cs_exps));

for k = 1:length(cs_exps)
   h = semilogx(ax1, freqs, 10*log10(abs(UZ_jump_s{k}))); 
   hold(ax1, 'on');
   h.DisplayName = name_s{k};
   grid(ax1, 'on')
   hands(k) = h;
   
   semilogx(ax2, freqs, 10*log10(abs(UZ_pre_jump_s{k})));
   hold(ax2, 'on');
   grid(ax2, 'on')
end
title(ax1, 'at jump')
title(ax2, 'pre-jump')
ylm = ylim(ax1);
ylim(ax2, ylm)
legend(hands);








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


function [UZ, UZ_pre_jump, UZ_jump, freqs] = psd_around_jumps(cs_exp, skips)
    N_mu = length(cs_exp.idx_state_s.scan);
    
    UZ = 0;
    UZ_pre_jump = 0;
    UZ_jump = 0;
    
    n_pre_jump=0;
    n_jump=0;

    for k = 1:N_mu
        idxk = cs_exp.idx_state_s.scan{k};
        uzk = cs_exp.uz(idxk);
        
        [UZK, freqs] = power_spectrum_local(uzk); %, AFM.Ts);
        
        UZ = UZK + UZ;
        
        if ~isempty(intersect(k, skips-1))
            UZ_pre_jump = UZK + UZ_pre_jump;
            n_pre_jump = n_pre_jump + 1;
        end
        
        if ~isempty(intersect(k, skips))
            UZ_jump = UZK + UZ_jump;
            n_jump = n_jump + 1;
        end
    end
    UZ = UZ/N_mu;
    UZ_pre_jump = UZ_pre_jump/n_pre_jump;
    UZ_jump = UZ_jump/n_jump;
end

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



