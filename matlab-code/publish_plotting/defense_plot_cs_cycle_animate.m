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
exp_date = '4-26-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-150prescan-notconnect_out_4-26-2019-01.csv',...
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
% taken from cycle 3, where tip detaches.
cs_exp1 = cs_exps{1};
zdfl_free = -0.637;
cs_exp1.ze = cs_exp1.ze - zdfl_free;


cs_exp1.x = cs_exp1.x*AFM.volts2mic_xy;
cs_exp1.y = cs_exp1.y*AFM.volts2mic_xy;
cs_exp1.uz = cs_exp1.uz*AFM.volts2nm_z;

%%
[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps{1}.plot_all_cycles(axs{1:4});



%%
clc
state = {'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13];
    'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up';
    'move', 'tdown', 'tsettle', 'scan', 'tup'};

Fig = mkfig(100, 4, 7); clf
ha = tight_subplot(3, 1, [.05, .05], [0.075, .03], [.12, .01], false);


tstart = 14.65;
CS_idx1 = cs_exps{1}.find_cycle_idx(tstart);
CS_idx2 = CS_idx1 + 1;



xlm1 = [14.91, 15.1];
hands = gobjects(length(state_names)+1, 3);

ha(1).XTickLabel = [];
ha(2).XTickLabel = [];



idxs = [4, 5, 1, 2, 3, 4];
cs_idxs = [CS_idx1, CS_idx1, CS_idx1+1, CS_idx1+1, CS_idx1+1, CS_idx1+1];
for k=1:6
hands(k, :) = plot_trajs(cs_exp1, state(:, idxs(k)), cs_idxs(k), ha);

% hands(2) = plot_trajs(cs_exp1, state(:, 5), CS_idx1, ha);
% hands(3) = plot_trajs(cs_exp1, state(:, 1), CS_idx2, ha);
% hands(4) = plot_trajs(cs_exp1, state(:, 2), CS_idx2, ha);
% hands(5) = plot_trajs(cs_exp1, state(:, 3), CS_idx2, ha);
end
% h_next = plot_trajs(cs_exp1, state(:, 4), CS_idx2, ha);

ylabel(ha(1), 'X [$\mu$m]', 'FontSize', 12)
ylabel(ha(2), '$u_Z$ [nm]', 'FontSize', 12)
ylabel(ha(3), 'deflection [v]', 'FontSize', 12)

ylim(ha(3), [-0.4, 0.45])


legend(hands(1:end-1,3), 'FontSize', 12, 'Location', 'SouthEast')

%%

% Why the fuck does matlab rescale my figure when I un-show a plot????
for j=1:3
    ylm = ylim(ha(j));
    ylim(ha(j), ylm);
    
    xlm = xlim(ha(j));
    xlim(ha(j), xlm);
end

%%

for k=2:6
    for j=1:3
        hands(k, j).Visible = 'Off';
    end
end


%     \includegraphics<1-2>[width=0.9\textwidth]{figures/CS_cycle_anime-01.pdf}
%     \includegraphics<3>[width=0.9\textwidth]{figures/CS_cycle_anime-02.pdf}
%     \includegraphics<4>[width=0.9\textwidth]{figures/CS_cycle_anime-03.pdf}
%     \includegraphics<5>[width=0.9\textwidth]{figures/CS_cycle_anime-04.pdf}
%     \includegraphics<6>[width=0.9\textwidth]{figures/CS_cycle_anime-05.pdf}

base_name = 'CS_cycle_anime_NEW-'; %01.pdf
for k=1:6
    fig_name = sprintf('%s0%d', base_name, k);
    for j=1:3
       hands(k, j).Visible = 'On'; 
    end
    
    save_fig(Fig, fullfile(PATHS.defense_fig(), fig_name), true)
end



%%


leg = legend(hands, 'Position', [0.0933 0.6311 0.1543 0.1000]);

% save_fig(Fig, fullfile(PATHS.thesis_root, 'plots-afm-cs-final/figures/cs_cycle'))


function h = plot_trajs(cs_exp1, state, CS_idx1, ha)
    traj_names = {'x', 'uz', 'ze'};
            
    %     names = {'$X$ [$\mu$m]', '$u_Z [nm]$', 'deflection [v]'};
    h = gobjects(1, 3);
    for j=1:length(traj_names)
        h_j = ha(j);
        
        h(j) = cs_exp1.plot_traj_from_csidx_by_state(CS_idx1, CS_idx1,...
            state{3}, traj_names{j}, h_j, 0, 'color', state{1});
        h(j).DisplayName = state{2};
        grid(h_j, 'on')

    end

end

