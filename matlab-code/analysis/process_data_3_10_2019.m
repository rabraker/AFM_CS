clc
clear
%
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
exp_date = '3-11-2019-KI-nosched'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-11-2019-01.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-11-2019-01.csv',...
'raster_scan_512pix_5mic_02Hz_out_3-11-2019-01.csv',...
'raster_scan_512pix_5mic_04Hz_out_3-11-2019-01.csv',...
'raster_scan_512pix_5mic_05Hz_out_3-11-2019-02.csv',...
'raster_scan_512pix_5mic_06Hz_out_3-11-2019-01.csv',...
};

cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-isconnect_out_3-11-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-250prescan-isconnect_out_3-11-2019-01.csv',...
'cs-traj-512pix-7perc-500nm-5mic-02Hz-250prescan-isconnect_out_3-11-2019-01.csv'  
};



for k=6:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
end

clc
x1s = [30, 32, 35, 37, 37, 86];
x2s = [499, 498, 498 452, 448, 449];
figbase = 10;
for k=1:length(rast_exps)
  rast_exp2 = copy(rast_exps{k});
%   uz = detrend(rast_exp2.uz);
  rast_exp2.bin_raster_really_slow(@detrend);
  
  pixmats_raw{k} = rast_exp2.pix_mat(1:end-1, 1:end-1);
%   pixmats{k} = pixmats_raw{k};
  pixmats{k} = pin_along_column(pixmats_raw{k}, x1s(k), x2s(k));
  pixmats{k} = pixmats{k} - mean(pixmats{k}(:));
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(pixmats{k}, figbase*k, stit)
%   stit = sprintf('(raw) scan %d', k);
%   plot_raster_data(pixmats_raw{k}, (figbase-5)*k, stit)
  
end

%%
data_root = PATHS.cs_image_data(size_dir, exp_date);

for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
%   cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth);
  cs_exps{k}.print_state_times();
end


[~, axs] = make_cs_traj_figs(figbase, 2);
cs_exps{1}.plot_all_cycles(axs{1:2});



bp = true;
for k=1:length(cs_exps)
  cs_exps{k}.process_cs_data(false, []);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
    cs_exps{k}.solve_bp(true, true);
    
    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, %% %.2f sampling\nTotal time: %.2f',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    f10=mkfig(1001 + 2*k, 6, 7.5); clf
    ax4 = subplot(3,1,[1,2]);
    ax4_2 =subplot(3,1,3);
    
    ImshowDataView.setup(f10);
    cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
    ImshowDataView.imshow(cs_exps{k}.Img_bp, [], ax4, ax4_2, cb_exp)
    title(ax4, stit)
  end

end
%%
mu = 75
Img_filts = {};
for k=1:length(cs_exps)
Img_filts{k} = cs_exps{k}.Img_bp - min(cs_exps{k}.Img_bp(:));

Img_filts{k} = SplitBregmanROF(Img_filts{k}, mu, 0.001);
figure(101+k)
mesh(Img_filts{k}), colormap('gray')
end

slice = 20:512-20;

im_master = pixmats{1} - mean(pixmats{1}(:));
% im_master = SplitBregmanROF(im_master, mu, 0.001);
fprintf('---------------------------------------------------\n');
scan_metrics = {};
scan_metrics{1} = ScanMetrics('ssim', 1, 'psnr', Inf,...
  'quality', rast_exps{1}.damage_metric(),...
  'damage', rast_exps{1}.quality_metric(),...
  'rate', rast_exps{1}.meta_in.raster_freq,...
  'time', rast_exps{1}.time_total,...
  'coverage', 100,...
  'type', 'raster');

for k=2:length(rast_exps)
  imk = pixmats{k} - mean(pixmats{k}(:));
  imk = SplitBregmanROF(imk, mu, 0.001);
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);
  damage = rast_exps{k}.damage_metric();
  quality = rast_exps{k}.quality_metric();
  
  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', rast_exps{k}.meta_in.raster_freq,...
    'time', rast_exps{k}.time_total,...
    'coverage', 100,...
    'type', 'raster');
  
  
  fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);

  scan_metrics{end+1} = csm;
end


fprintf('---------------------------------------------------\n');
for k=1:length(cs_exps)
  imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
  imk = SplitBregmanROF(imk, mu, 0.001);
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);

  damage = cs_exps{k}.damage_metric();
  quality = cs_exps{k}.quality_metric();

  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
    'time', length(cs_exps{k}.x)*cs_exps{k}.Ts,...
    'coverage', cs_exps{k}.meta_in.actual_sub_samble_perc,...
    'type', 'CS');
  
  fprintf('(CS     %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);
  
 
  scan_metrics{end+1} = csm;
end
%%
S = scan_metrics_table(scan_metrics)
fid = fopen('notes/cs_raster_table_3-11-2019_01.tex', 'w+');
fprintf(fid, '%s', S);
fclose(fid);

S = state_times_table(cs_exps)


fid = fopen('notes/cs_state_times_table_3-11-2019_01.tex', 'w+');
fprintf(fid, '%s', S);
fclose(fid);


function S = scan_metrics_table(csm_s)
  S = 'type &  rate (Hz) & PSNR & SSIM & time [s] & damage\\';
  S = sprintf('%s\n\\toprule\n', S);
  for k=1:length(csm_s)
    csm = csm_s{k};
    S = sprintf('%s%s & %.2f & %.2f & %.2f & %.1f & %.2f\\\\\n', S,...
      csm.type, csm.rate, csm.psnr, csm.ssim, csm.time, csm.damage);
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
