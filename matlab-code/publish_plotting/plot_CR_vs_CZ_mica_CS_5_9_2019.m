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

data_root = '/media/labserver/afm-cs/z-bounce/5-9-2019/';

%   '/media/labserver/afm-cs/imaging/cs-imaging/5microns/3-28-2019/'}
% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-9-2019-01.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-9-2019-02.csv',...
};

cs_exps = {};
for k=1:2 %length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
%   gg = zpk([], [], 1, AFM.Ts);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

for k=1:length(cs_exps)
   cs_exps{k}.uz  = cs_exps{k}.uz*AFM.volts2nm_z;
end

%%
[~, axs1] = make_cs_traj_figs(figbase, 2);
cs_exps{1}.plot_all_cycles(axs1{1:2});

[~, axs2] = make_cs_traj_figs(figbase+100, 4);
cs_exps{2}.plot_all_cycles(axs2{1:4}, [], 512);

linkaxes([axs1{:}, axs2{:}], 'x')
%%

skips = [9:8:160*5];
ii_cz=1;
ii_cr=2;
[UZ_cz, UZ_pre_jump_cz, UZ_jump_cz, freqs_cz] = psd_around_jumps(cs_exps{ii_cz}, skips);
[UZ_cr, UZ_pre_jump_cr, UZ_jump_cr, freqs_cr] = psd_around_jumps(cs_exps{ii_cr}, skips);

width = 3.5;
height = 3;
gap = 0.1;
margh = [0.13, .025];
margw = [0.15, .025];


% ---------------------------------- %
F2 = mkfig(20, width, height); clf
ax2 = tight_subplot(1, 1, gap, margh, margw, false);

h0 = semilogx(ax2, freqs, 10*log10(abs(UZ_pre_jump_cz)));
h0.DisplayName = 'choose-$\zeta$'
hold(ax2, 'on')
grid(ax2, 'on')
h1 = semilogx(ax2, freqs, 10*log10(abs(UZ_pre_jump_cr)));
h1.DisplayName = 'constant-$\rho$'
grid(ax2, 'on')
xlabel(ax2, 'Hz')
ylabel(ax2, 'PSD')
leg1 = legend([h0, h1], 'Position', [0.1492 0.1286 0.5084 0.1712], 'FontSize', 14);

% ---------------------------------- %
F23 = mkfig(23, width, height); clf
ax3 = tight_subplot(1, 1, gap, margh, margw, false);

h2 = semilogx(ax3, freqs, 10*log10(abs(UZ_jump_cz)));
hold(ax3, 'on')
h2.DisplayName = 'choose-$\zeta$'

h3 = semilogx(ax3, freqs_b, 10*log10(abs(UZ_jump_cr)));
h3.DisplayName = 'constant-$\rho$'



leg2 = legend([h2, h3], 'Position', [0.1481 0.1288 0.7262 0.1712], 'FontSize', 14);
grid(ax3, 'on')
xlabel(ax3, 'Hz')
ylabel(ax3, 'PSD')




%%



F25 = mkfig(25, width, height*1.5); clf
axs2 = tight_subplot(2, 1, gap, [0.1, 0.1], margw, false);

% title(axs1(1), '$xy$ fixed', 'FontSize', 14)
tstart_cr = 1.45; 
tend_cr = 1.6; 
zbounce = false;
plot_uz_time_range(cs_exps{ii_cr}, axs2(1), tstart_cr, tend_cr, zbounce, 'uz');
zbounce = false;
h = plot_uz_time_range(cs_exps{ii_cr}, axs2(2), tstart_cr, tend_cr, zbounce, 'x');


grid(axs2(1), 'on')
grid(axs2(2), 'on')

% ylim(axs2(1), [-5, 20])


xlim(axs2(1), [tstart_cr, tend_cr])
xlim(axs2(2), [tstart_cr, tend_cr])

ylim(axs2(2), [-0.5, 5])
ylim(axs2(1), [-5, 100])

ylabel(axs2(1), '$u_Z$~[nm]', 'FontSize', 14)
ylabel(axs2(2), '$x$~[$\mu$m]', 'FontSize', 14)
xlabel(axs2(2), 'time [s]', 'FontSize', 14)

legend(h, 'location', 'northeast');

% ----------------------- choose-zeta -----------------
F26 = mkfig(26, width, height*1.5); clf
axs6 = tight_subplot(2, 1, gap, [0.1, 0.1], margw, false);

title(axs2(1), 'constant-$\rho$ (BW = 111~Hz)', 'FontSize', 14)
tstart_cz = 1.49; 
tend_cz = 1.65; 

zbounce = false;
plot_uz_time_range(cs_exps{ii_cz}, axs6(1), tstart_cz, tend_cz, zbounce, 'uz');
zbounce = false;
h = plot_uz_time_range(cs_exps{ii_cz}, axs6(2), tstart_cz, tend_cz, zbounce, 'x');


grid(axs6(1), 'on')
grid(axs6(2), 'on')
title(axs6(1), 'choose-$\zeta$ (BW = 62~Hz)', 'FontSize', 14)

xlim(axs6(1), [tstart_cz, tend_cz])
xlim(axs6(2), [tstart_cz, tend_cz])

ylim(axs6(2), [-0.5, 5])
ylim(axs6(1), [-5, 100])
ylabel(axs6(1), '$u_Z$~[nm]', 'FontSize', 14)
ylabel(axs6(2), '$x$~[$\mu$m]', 'FontSize', 14)
xlabel(axs6(2), 'time [s]', 'FontSize', 14)

legend(h, 'location', 'northeast');
%%
save_fig(F23, fullfile(PATHS.defense_fig(), 'cr_vs_cz_mica_jump_PSD'), true)
save_fig(F25, fullfile(PATHS.defense_fig(), 'cr_mica_jump_xuz'), true)
save_fig(F26, fullfile(PATHS.defense_fig(), 'cz_mica_jump_xuz'), true)

%%

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



