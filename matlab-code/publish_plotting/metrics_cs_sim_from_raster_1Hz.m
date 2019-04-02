clc
clear

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
exp_date = '3-21-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_01Hz_out_3-21-2019-01.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-21-2019-03.csv',...
};



rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
end


x1s = [44,   49];
x2s = [459, 459];
figbase = 10;
for k=1:length(rast_exps)
  rast_exps{k}.bin_raster_really_slow(@detrend);
  
  pixmats_raw{k} = rast_exps{k}.pix_mat(1:end, 1:end);
  pixmat_ = pin_along_column(rast_exps{k}.pix_mat, x1s(k), x2s(k));
  
  rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:));
  rast_exps{k}.pin_idx_s = [x1s(k), x2s(k)];
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit)
end

for k=1:length(rast_exps)
  rast_exps{k}.save()
end
%%

% Cs simulations using raster data.
master_idx = 1;
rfreq_master = rast_exps{master_idx}.meta_in.raster_freq;
im_rast_sim = rast_exps{master_idx}.pix_mat_pinned;
im_rast_sim = im_rast_sim - mean(im_rast_sim(:));
use_dct2=false;

[n, m] = size(im_rast_sim);
opts = l1qc_dct_opts('l1_tol', 0.0001, 'verbose', 1);
cs_idx = [1, 5];
rast_idx = [2,3];
rng(1);
width =5;    % microns
npix = 512;  % image resolution.
mu_length = 0.5;  % 500 nm. length of the horizontal mu-path.
pix_per_micron = npix/width;
mu_pix = ceil(mu_length*pix_per_micron);

sub_sample_frac_s = [0.075, 0.10, 0.15, 0.2, 0.25];

cs_sims = {};


for k=1:length(sub_sample_frac_s)
  
  sub_sample_frac = sub_sample_frac_s(k)
  [pix_mask] = mu_path_mask(mu_pix, npix, npix, sub_sample_frac, false);

  figure(4)
  imagesc(pix_mask)
  
  cs_sims{k} = CsSim(im_rast_sim, pix_mask);
  fprintf('Solving BP problem for %.2f sampling...', sub_sample_frac*100);
  cs_sims{k}.solve_bp(false, use_dct2, opts);
  fprintf('done\n');
end

mu = Inf;
F1 = mkfig(3000, 7, 4.55); clf
[ha, pos] = tight_subplot(2, 3, [.05, .01 ], [.01, .04], .015);
thresh = (20/7)*(1/1000)*20;

imagesc(ha(1), im_rast_sim, [-thresh, thresh])
title(ha(1), sprintf('master (%.1f Hz)', rfreq_master))
set(ha(1), 'YTick', [], 'XTick', [])

for k=1:length(sub_sample_frac_s)
  imk = cs_sims{k}.Img_bp;
  if ~isinf(mu)
    imk = SplitBregmanROF(imk, mu, 0.001);
  end
%   im1_ontok_fit = norm_align(, im_rast_sim);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im_rast_sim, cs_sims{k}.Img_bp);
  
  stit = sprintf('\\%% %.2f: psnr=%.2f, ssm=%.2f',...
    cs_sims{k}.sample_frac()*100, psn_1k, ssm_1k);
  fprintf('%s\n', stit);
%   cs_exps{cs_idx(k)}.meta_in.actual_sub_samble_perc
 
  imagesc(ha(k+1), imk, [-thresh, thresh])
  colormap('gray')
  title(ha(k+1), stit)
  set(ha(k+1), 'YTick', [], 'XTick', [])
end
%%
save_fig(F1, fullfile(PATHS.thesis_root, 'plots-afm-cs-final/figures/cs_sim_1Hz_raster'), false)



%% ---------------------------------------------------------------------------
im_rast_sim2 = rast_exps{2}.pix_mat_pinned;
im_rast_sim2 = im_rast_sim2 - mean(im_rast_sim2(:));

for k=1:length(sub_sample_frac_s)
  
  sub_sample_frac = sub_sample_frac_s(k)
  [pix_mask] = mu_path_mask(mu_pix, npix, npix, sub_sample_frac, false);

  figure(4)
  imagesc(pix_mask)
  
  cs_sims2{k} = CsSim(im_rast_sim2, pix_mask);
  fprintf('Solving BP problem for %.2f sampling...', sub_sample_frac*100);
  cs_sims2{k}.solve_bp(false, use_dct2, opts);
  fprintf('done\n');
end
%
mu = Inf;
slice = 25:1:cs_sims2{k}.npix-25;
F2 = mkfig(3001, 7, 4.5); clf
[ha, pos] = tight_subplot(2, 3, [.05, .01 ], [.01, .04], .015);
thresh = (20/7)*(1/1000)*20;

imagesc(ha(1), im_rast_sim, [-thresh, thresh])
title(ha(1), sprintf('master (%.1f Hz)', rfreq_master))
set(ha(1), 'YTick', [], 'XTick', [])

for k=1:length(sub_sample_frac_s)
  imk = cs_sims2{k}.Img_bp;
  if ~isinf(mu)
    imk = SplitBregmanROF(imk, mu, 0.001);
  end
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_rast_sim);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);
  
  stit = sprintf('\\%% %.2f: PSNR=%.2f, SSIM=%.2f',...
    cs_sims2{k}.sample_frac()*100, psn_1k, ssm_1k);
  fprintf('%s\n', stit);
%   cs_exps{cs_idx(k)}.meta_in.actual_sub_samble_perc
 
  imagesc(ha(k+1), imk, [-thresh, thresh])
  colormap('gray')
  title(ha(k+1), stit)
  set(ha(k+1), 'YTick', [], 'XTick', [])
  
end
%%
imk = im_rast_sim2;
imk_slice = imk(slice, slice);
im1_ontok_fit = norm_align(imk_slice, im_rast_sim);
[psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);
  
fprintf('(master to 1hz) PSNR=%.2f, SSIM=%.2f\n', psn_1k, ssm_1k);

%%
save_fig(F2, fullfile(PATHS.thesis_root,'plots-afm-cs-final/figures/cs_sim_1Hz_raster_from1Hz'), false)



% ---------------------------------------------------------------------------- %
% ---------------------- Local Functions ------------------------------------- %
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