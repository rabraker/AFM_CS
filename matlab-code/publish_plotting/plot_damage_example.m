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
%

cs_exp.meta_exp.z_axis_params.enable_KI_sched

zdfl_free = -0.5479;
cs_exp.ze = cs_exp.ze - zdfl_free;
%%
figbase = 20;
figure(1);
ax1 = gca();
Fig = mkfig(figbase+2, 5, 3); clf
[ha, pos] = tight_subplot(1, 1, [.01], [.15, .02], [.12, .01], false)

hands = cs_exp.plot_all_cycles(ax1, ha);

% cs_exp.plot_traj_in_time_interval_by_state(tstart, tend, 'ze'

ref = cs_exp.meta_exp.z_axis_params.setpoint_scan;


ylm = [-0.9478    0.4434] - zdfl_free;
set(ha, 'XLim', xlm, 'YLim', ylm)

tstart = 3.8566;
tend = 4.2150;
xlm = [tstart, tend];


idx_set = cs_exp.get_idx_by_state_in_time_range('tsettle', tstart, tend);

% idx_scan = cs_exp.get_idx_by_state_in_time_range('scan', tstart, tend);
idx_scans = get_idx_scans_in_time_range(cs_exp, tstart, tend+.2)

for k=1:length(idx_set)
    idx_k = idx_set{k};
    tt = cs_exp.t(idx_k);
    zz = cs_exp.ze(idx_k);
    
    idx_pos = find(zz > ref_eff);
    h_no = ciplot(tt(idx_pos)*0+ref_eff, zz(idx_pos), tt(idx_pos), 'r', ha);
    alpha(h_no, '0.25')
end

for k=1:length(idx_scans)
    idx_k = idx_scans{k};
    tt = cs_exp.t(idx_k);
    zz = cs_exp.ze(idx_k);
    
    idx_pos = find(zz > ref_eff);
    h_no = ciplot(tt(idx_pos)*0+ref_eff, zz(idx_pos), tt(idx_pos), 'r', ha);
    alpha(h_no, '0.25')
end

% h_no = ciplot([ref, ref]-zdfl_free, [ylm(2), ylm(2)], xlm, 'r', ha);
% h_ok.DisplayName = 'No Penalty';
h_no.DisplayName = 'Penalize';
xlabel(ha, 'time [s]')
ylabel(ha, 'deflection [V]')
title(ha, '')

leg = legend(hands(:, 3));
set(leg, 'NumColumns', 2, 'location', 'northwest')

xlm = [3.96, 4.15];
xlim(xlm)


ha.XLabel.FontSize = 14;
ha.YLabel.FontSize = 14;
ha.Legend.FontSize = 12;
set(leg, 'NumColumns', 1);
leg.Position = [0.7097 0.5648 0.2784 0.4226];
%%
save_fig(Fig, fullfile(PATHS.thesis_fig_final, 'damage_illustration'))
save_fig(Fig, fullfile(PATHS.defense_fig(), 'damage_illustration'), true)
%%

exp_date = '4-26-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_01Hz_out_4-26-2019-01.csv',...
'raster_scan_512pix_5mic_10Hz_out_4-26-2019-01.csv',...
};


%%
rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths, 'load_full', true);
  if rast_exps{k}.time_total == 0
      rast_exps{k}.time_total = rast_exps{k}.samps_per_period*rast_exps{k}.npix*AFM.Ts;
  end
end


%%


Fig = mkfig(figbase+2, 7, 3); clf
[ha, pos] = tight_subplot(1, 1, [.01], [.15, .08], [.12, .02], false);




row_idx = 174;

rates = [1, 10];
hands = gobjects(2,1);
names = {'1 Hz', '10 Hz'};
for k=1:2
    
    idx = rast_exps{k}.samps_per_period*row_idx : (rast_exps{k}.samps_per_period*(row_idx+1) - rast_exps{k}.samps_per_line);
    t = (0:length(rast_exps{k}.ze)-1')*AFM.Ts;
        zz = rast_exps{k}.ze(idx);
    idx_pos = find(zz > -0.3);
    tt = t(idx)*rates(k) - t(idx(1))*rates(k);
    
    hands(k) = plot(ha, tt, zz);
    hands(k).DisplayName = names{k};
    

    h_no = ciplot(tt(idx_pos)*0-0.3, zz(idx_pos), tt(idx_pos), 'r');
    alpha(h_no, '0.25')
    
    hold(ha, 'on')
    ylim(ha, [-0.45, -0.15])
    xlim(ha, [0, (+0.5)])
end

leg = legend(hands, 'Position', [0.7906 0.7725 0.1905 0.1490], 'FontSize', 12);

ylabel('deflection [v]', 'FontSize', 14)
xlabel('normalized time', 'FontSize', 14)
title('raster', 'FontSize', 14)
save_fig(Fig, fullfile(PATHS.defense_fig(), 'damage_raster_illustration'), true)


function idx_scans = get_idx_scans_in_time_range(cs_exp, tstart, tend)
    
    scan_idx_cell = cs_exp.idx_state_s.scan;
    idx_scans = {};
   for k=1:length(scan_idx_cell) 
    
       t_kk = cs_exp.t(scan_idx_cell{k});
       if tstart < t_kk(1) && tend > t_kk(end) 
           idx_scans{end+1} = scan_idx_cell{k};
       end
       if all(cs_exp.t(scan_idx_cell{k}) > tend)
           break;
       end
   end
end







