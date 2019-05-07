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

data_root_s = '/media/labserver/afm-cs/z-bounce/4-19-2019/';

%   '/media/labserver/afm-cs/imaging/cs-imaging/5microns/3-28-2019/'}
% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
    'cs-traj-512pix-3perc-500nm-5mic-01Hz-50prescan-notconnect_out_4-19-2019-02.csv',...
    'cs-traj-512pix-3perc-500nm-5mic-01Hz-50prescan-notconnect_out_4-19-2019-03.csv',...
    'cs-traj-512pix-3perc-500nm-5mic-01Hz-50prescan-notconnect_out_4-19-2019-04.csv'
};

% cs_files = {...
% 'cs-traj-z-bounce_1500samples_Npath160_prescan250_out_4-18-2019-03.csv',...
% 'cs-traj-z-bounce_1500samples_Npath160_prescan250_out_4-18-2019-04.csv',...
% };

cs_exps = {};
for k=1:3 %length(cs_files)
  cs_paths = get_cs_paths(data_root_s, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
%   gg = zpk([], [], 1, AFM.Ts);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  cs_exps{k}.uz = cs_exps{k}.uz - min(cs_exps{k}.uz);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end
%%

%
[~, axs1] = make_cs_traj_figs(figbase, 2);
cs_exps{1}.plot_all_cycles(axs1{1:2});

[~, axs3] = make_cs_traj_figs(figbase+200, 3);
cs_exps{2}.plot_all_cycles(axs3{1:3});

% [~, axs2] = make_cs_traj_figs(figbase+100, 4);
% cs_exps{2}.plot_all_cycles(axs2{1:4}, [], 512);

% [~, axs2] = make_cs_traj_figs(figbase+100, 4);
% cs_exps{2}.plot_all_cycles(axs2{1:4}, [], 512);


linkaxes([axs1{:}, axs3{:}], 'x')
%%


width = 6;
height = 6;
gap = 0.05;
margh = 0.1;
margw = 0.1;
margh = [0.1, .05];
margw = [0.08, .04];

F1 = mkfig(21, width, height); clf
ax1 = tight_subplot(2, 2, gap, margh, margw, false);
ax1 = reshape(ax1', 2, 2)';
%
% F2 = mkfig(22, width, height); clf
% ax2 = tight_subplot(2, 1, gap, margh, margw, false);

%
to = 1;
tstart1 = 2.12; %3.46 + to;
tend1 = 2.24; %3.5 + to;

zbounce = false;
ze_sub = -0.448;
pre_ext = 150;
plot_uz_time_range(cs_exps{1}, ax1(1,1), tstart1, tend1, zbounce, 'uz', 4/AFM.volts2nm_z, pre_ext);
h = plot_uz_time_range(cs_exps{1}, ax1(2, 1), tstart1, tend1, zbounce, 'ze', ze_sub, pre_ext)
%
% plot_uz_time_range(cs_exps{2}, ax1(1, 2), tstart1, tend1, zbounce, 'uz');
% plot_uz_time_range(cs_exps{2}, ax1(2, 2), tstart1, tend1, zbounce, 'ze');

plot_uz_time_range(cs_exps{3}, ax1(1, 2), tstart1, tend1, zbounce, 'uz', 0, pre_ext);
plot_uz_time_range(cs_exps{3}, ax1(2, 2), tstart1, tend1, zbounce, 'ze', ze_sub, pre_ext);

%
for j=1:2
    for k=1:2
        grid(ax1(k, j), 'on')
        
    end
    xlabel(ax1(2,j), 'time [s]', 'FontSize', 14)
end
ylim(ax1(1, :), [-0.06, 1.2]*AFM.volts2nm_z)
% ylim(ax1(1), [-0.06, 1.2]*AFM.volts2nm_z)

title(ax1(1, 1), 'big lift', 'FontSize', 14)
% title(ax1(1, 2), "medium lift", 'FontSize', 14)
title(ax1(1, 2), "small lift ", 'FontSize', 14)

% delt = 
delt = 0.081;
tstart1 = 2.12;
tstart2 = 2.122;
tstart3 = 2.122;

xlim(ax1(:,1), [tstart1, tstart1+delt])
% xlim(ax1(:,2), [tstart2, tstart3+delt])
xlim(ax1(:,2), [tstart3, tstart3+delt])




ylabel(ax1(1), '$u_z$ [nm]', 'FontSize', 14)
ylabel(ax1(2), '$Z_{\textrm{dfl}}$ [v]', 'FontSize', 14)
% ylabel(ax2(1), '$u_z$  [nm]', 'FontSize', 14)
% ylabel(ax2(2), '$Z_{\textrm{dfl}}$ [v]', 'FontSize', 14)


legend(h, 'location', 'southwest', 'FontSize', 12)


% save_fig(F2, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/prescan_uz_example'), false)
%%
if 1
    save_fig(F1, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/justify_small_lift'), false)
    save_fig(F1, fullfile(PATHS.defense_fig(), 'justify_small_lift'), true)
end
ylim(ax1(1, :), [-0.0, 0.1]*AFM.volts2nm_z)


% xlim(ax1, [2.1, 2.28])
% xlim(ax2, [2.1, 2.28])

if 1
    save_fig(F1, fullfile(PATHS.thesis_root(), 'plots-afm-cs-final/figures/justify_small_lift_zoom'), false)
    save_fig(F1, fullfile(PATHS.defense_fig(), 'justify_small_lift_zoom'), true)
end
%%
exp_date = '3-20-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-20-2019-01.csv',...
};


data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps2 = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps2{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  cs_exps2{k}.print_state_times();
  cs_exps2{k}.sub_sample_frac()
end
%%
Fig = mkfig(100, width, width); clf
ha = tight_subplot(2, 1, [.05, .05], [0.075, .03], [.1, .01], false);
CS_idx1 = cs_exps{1}.find_cycle_idx(16.31);
CS_idx2 = CS_idx1 + 3;
tstart4 = 16.45
ze_sub = -0.637;
delt2 = 0.38;
plot_uz_time_range(cs_exps2{1}, ha(1), tstart4, tstart4+delt2, false, 'uz', 0, 0)
hands = plot_uz_time_range(cs_exps2{1}, ha(2), tstart4, tstart4+delt2, false, 'ze', ze_sub, 0);

legend(hands, 'location', 'southeast', 'FontSize', 12)

ylabel(ha(1), '$U_Z$ [nm]', 'FontSize', 14)
ylabel(ha(2), 'deflection [v]', 'FontSize', 14)
xlabel(ha(2), 'time [s]', 'FontSize', 14)

xlim(ha, [16.45, 16.35+delt2])
grid(ha(1), 'on')
grid(ha(2), 'on')

%%
save_fig(Fig, fullfile(PATHS.defense_fig(), 'small_lift_justify_ze_norm'), true)

%%
function hands =  plot_uz_time_range(cse, ax, tstart, tend, zbounce, signal, ze_sub, pre_ext)
indc = {'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
        'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};

    hands = gobjects(5, 1);
    hold(ax, 'on')
    state_s = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
    if signal == 'uz'
        scl = AFM.volts2nm_z;
    elseif signal == 'x'
        scl = AFM.volts2mic_xy;
    else
        scl = 1;
    end
    
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
            uzk = cse.(signal)(idxs{j}) - ze_sub;
            
            if k ~=  4 && (pre_ext >= 1)
                hands(k) =    plot(ax, t, uzk*scl, 'Color', clr);
                hands(k).DisplayName = indc{2, k};
            else
               hands(k) = plot(ax, t(pre_ext+1:end), uzk(pre_ext+1:end)*scl, 'Color', clr);
               hands(k).DisplayName = indc{2, k};
                plot(ax, t(1:pre_ext), uzk(1:pre_ext)*scl, 'Color', indc{1, 3});
            end
        end
        
        
    end
end


function [Pxx, freqs] = power_spectrum_local(X)

  N = length(X);

  [Pxx, freqs] = periodogram(detrend(X), [], N, 1/AFM.Ts);
end



