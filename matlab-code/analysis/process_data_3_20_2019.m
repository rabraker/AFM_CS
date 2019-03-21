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
raster_files = {...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_04Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_05Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_08Hz_out_3-20-2019-01.csv',...
'raster_scan_512pix_5mic_10Hz_out_3-20-2019-01.csv',...
};
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


%%
rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
end
%%

x1s = [30, 47, 53, 51, 59, 51];
x2s = [496, 465, 464, 468, 472, 465];
figbase = 10;
for k=1:length(rast_exps)
  rast_exps{k}.bin_raster_really_slow(@detrend);
  
  pixmats_raw{k} = rast_exps{k}.pix_mat(1:end, 1:end);
%   rast_exps{k}.pix_mat_pinned = pixmats_raw{k};
  pixmat_ = pin_along_column(rast_exps{k}.pix_mat, x1s(k), x2s(k));
  
  rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:));
  rast_exps{k}.pin_idx_s = [x1s(k), x2s(k)];
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit)
  % stit = sprintf('(raw) scan %d', k);
  % plot_raster_data(pixmats_raw{k}, (figbase-5)*k, stit)
end
%%
for k=1:length(rast_exps)
  rast_exps{k}.save()
end
%%
data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  % cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps{end}.plot_all_cycles(axs{1:4});


%%
bp = true;
use_dct2 = true;
thresh = (20/7)*(1/1000)*20;
opts = l1qc_dct_opts('l1_tol', 0.001);
for k=1:length(cs_exps)
  cs_exps{k}.process_cs_data(false, []);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
    cs_exps{k}.solve_bp(true, use_dct2, opts);

    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, \\%% %.2f sampling\nTotal time: %.2f',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{k}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    f10=mkfig(1001 + 2*k, 6, 7.5); clf
    ax4 = subplot(3,1,[1,2]);
    ax4_2 =subplot(3,1,3);
    
    ImshowDataView.setup(f10);
    cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
    ImshowDataView.imshow(cs_exps{k}.Img_bp-mean(cs_exps{k}.Img_bp(:)),...
      [-thresh, thresh], ax4, ax4_2, cb_exp)
    title(ax4, stit)
  end

  cs_exps{k}.save();
end
%%
mu = Inf;
Img_filts = {};
mxs = [];
thresh = (20/7)*(1/1000)*20;
DRng = 2*thresh;

for k=1:length(cs_exps)
Img_filts{k} = cs_exps{k}.Img_bp - min(cs_exps{k}.Img_bp(:));
if ~isinf(mu)
  Img_filts{k} = SplitBregmanROFAn(Img_filts{k}, mu, 0.001);
end
figure(101+k)
mx = max(Img_filts{k}(:)) - min(Img_filts{k}(:));
mxs = [mxs; mx];
mesh(Img_filts{k}), colormap('gray')
end

slice = 30:512-30;
master_idx = 2;
im_master = rast_exps{master_idx}.pix_mat_pinned - mean(rast_exps{master_idx}.pix_mat_pinned(:));
if ~isinf(mu)
  %im_master = SplitBregmanROF(im_master, mu, 0.001);
end
fprintf('---------------------------------------------------\n');
scan_metrics = {};
scan_metrics{1} = ScanMetrics('ssim', 1, 'psnr', Inf,...
  'quality', rast_exps{1}.damage_metric(),...
  'damage', rast_exps{1}.quality_metric(),...
  'rate', rast_exps{1}.meta_in.raster_freq,...
  'time', rast_exps{1}.time_total,...
  'coverage', 100,...
  'type', 'raster');

figure(3000);clf;
h = subplott(3,4);
fprintf('\n')
j=1;
for k=1:length(rast_exps)
  if k == master_idx
    continue;
  end
  imk = rast_exps{k}.pix_mat_pinned - mean(rast_exps{k}.pix_mat_pinned(:));
  if ~isinf(mu)
    %imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
%   [im1_ontok_fit, imk_slice] = align_by_metric(im_master, imk, [], 'psnr');
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);
  damage = rast_exps{k}.damage_metric();
  quality = rast_exps{k}.quality_metric();
  
  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', rast_exps{k}.meta_in.raster_freq,...
    'time', rast_exps{k}.time_total,...
    'coverage', 100,...
    'type', 'raster');
  
%   subplot(2,4,k-1)
  imshowpair(im1_ontok_fit, imk_slice, 'parent', h(j));
  title(h(j), sprintf('raster (%.1f Hz)', csm.rate));
  fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);

  scan_metrics{end+1} = csm;
  j=j+1;
end


fprintf('---------------------------------------------------\n');
cs_stats = {};
for k=1:length(cs_exps)
  imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);

  damage = cs_exps{k}.damage_metric();
  quality = cs_exps{k}.quality_metric();

  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
    'time', length(cs_exps{k}.x)*cs_exps{k}.Ts,...
    'coverage', cs_exps{k}.meta_in.actual_sub_samble_perc,...
    'type', 'CS');
  
  fprintf('(CS %.1f Hz, %% %.1f, %.1f) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    csm.rate, csm.coverage, cs_exps{k}.sub_sample_frac*100, psn_1k, ssm_1k, damage, csm.time);
  
%   figure(3000+k+length(rast_exps));
% subplot(2,4,k+length(rast_exps)-1)
  imshowpair(im1_ontok_fit, imk_slice, 'parent', h(k+length(rast_exps)-1));
  title(h(length(rast_exps)+k-1), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));
 
  scan_metrics{end+1} = csm;
  [t_cycle_avg, t_connect] = cs_exps{k}.estimate_mpt_connect_savings();
  cs_stats{k, 1} = csm;
  cs_stats{k,2} = t_cycle_avg;
  cs_stats{k,3} = t_connect;
end

%%




S1 = scan_metrics_table(scan_metrics)

S2 = state_times_table(cs_exps)
S3 = mpt_connect_table(cs_stats)

%%
fid = fopen('notes/tables/cs_raster_table_3-18-2019_mu100_dct2.tex', 'w+');
fprintf(fid, '%s', S1);
fclose(fid);

fid = fopen('notes/tables/cs_state_times_table_3-18-2019_mu100_dct2.tex', 'w+');
fprintf(fid, '%s', S2);
fclose(fid);




fid = fopen('notes/tables/cs_connect_table_3-18-2019_mu100_dct2.tex', 'w+');
fprintf(fid, '%s', S3);
fclose(fid);

function S = mpt_connect_table(cs_stats)
  S = 'coverage &  rate (Hz) & skipped cycle time (est) & connect time \\';
  S = sprintf('%s\n\\toprule\n', S);
  for k=1:length(cs_stats)
    csm = cs_stats{k, 1};
    t_cycle_est = cs_stats{k, 2};
    t_connect = cs_stats{k,3};
    S = sprintf('%s %.2f & %.2f & %.2f & %.1f \\\\\n', S,...
      csm.coverage, csm.rate, t_cycle_est, t_connect);
  end
end



function S = scan_metrics_table(csm_s)
  S = 'type &  rate (Hz) & PSNR & SSIM & time [s] & damage\\';
  S = sprintf('%s\n\\toprule\n', S);
  
  for k=1:length(csm_s)
    csm = csm_s{k};
    if strcmp(csm.type, 'CS')
      description = sprintf('CS (%.2f)\\%%', csm.coverage);
    else
      description = 'raster';
    end
    S = sprintf('%s%s & %.2f & %.2f & %.2f & %.1f & %.2f\\\\\n', S,...
      description, csm.rate, csm.psnr, csm.ssim, csm.time, csm.damage);
  end
end

function S = state_times_table(cs_exps)
  S = 'description & move & engage & tip-settle & scan & connect & tip-up & total\\';
  S = sprintf('%s\n\\toprule\n', S);
  for k=1:length(cs_exps)
    cst = cs_exps{k}.get_state_times();
    description = sprintf('%.1f Hz, %.2f\\%%',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc);
  
    S = sprintf('%s%s &%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\\n', S,...
      description, cst.move, cst.tdown, cst.tsettle, cst.scan,...
      cst.connect, cst.tup, cst.total);
  end
end




%%
function plot_raster_data(pixmat2, figbase, stit, plot_mesh)
  if nargin < 4
    plot_mesh = false;
  end
  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
%   F10 = figure(figbase+2); clf
  F10 = mkfig(figbase+2, 6, 7.5, false);
  ax1 = subplot(3,1,[1,2]);
  ax2 =subplot(3,1,3);
    

  lo = min(min(pixmat2));
  hi = max(max(pixmat2));
  
  
  imshow_dataview(flipud(pixmat2 - mean(pixmat2(:))), [-thresh, thresh], ax1, ax2)
  try
    grid(ax1, 'on')
  catch
    keyboard
  end
  
  colormap(ax1, 'gray')
  grid(ax2, 'on')
  ax1.GridAlpha = 1;
  ax2.GridAlpha = 1;
  title(ax1, stit)
  title(ax2, stit)
  if plot_mesh
    figure(figbase+4)
    mesh(pixmat2)
    xlabel('x')
    ylabel('y')
    title(stit)
  end
end
