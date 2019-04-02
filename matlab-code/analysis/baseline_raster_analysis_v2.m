% This file compares 8 slow raster scans all taken back to back at 0.5 Hz overa
% 5 micron area. The goal is to establish a baseline of what we should expect
% for the performance metrics.

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
exp_date = '3-13-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-09.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-10.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-11.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-12.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-13.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-14.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-15.csv',...
'raster_scan_512pix_5mic_5.00e-01Hz_out_3-13-2019-16.csv',...
};



pixmaps = cell(length(raster_files), 1);
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths);
end


npix = rast_exps{1}.npix;
Ts = rast_exps{1}.Ts;
width = rast_exps{1}.width;
stit = sprintf('Scan %d', k)

clc
% N_damage = 100;
% for k=1:length(rast_exps)
%   ze{k} = rast_exps{k}.ze;
%   ze{k} = resample(ze{k}, 1, 2); % half the sample rate
%   ze{k} = ze{k} - (-0.3); %mean(ze1);
%   ze_sort{k} = sort(ze{k}(ze{k}>0), 'descend');
%   %   sum(ze_sort{k}(1:N_damage))/N_damage
%   sum(ze{k}(ze{k}>0)) / length(ze{k}(ze{k}>0))
% end


%%
clc

x1s = [45,   45, 50,   50, 50,   52,  57, 55];
x2s = [450, 455, 455, 463, 419, 467. 474, 474];
figbase = 10;
for k=1:length(rast_exps)
  rast_exps{k}.bin_raster_really_slow(@detrend);
  
  pixmats_raw{k} = rast_exps{k}.pix_mat;
%   pixmats{k} = pixmats_raw{k};
  pixmats{k} = pin_along_column(pixmats_raw{k}, x1s(k), x2s(k));
  pixmats{k} = pixmats{k} - mean(pixmats{k}(:));
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(pixmats{k}, figbase*k, stit)
%   stit = sprintf('(raw) scan %d', k);
% if k==5
%   keyboard
% end
end
%%

for k=1:length(rast_exps)
  
  damage = rast_exps{k}.damage_metric();
  quality = rast_exps{k}.quality_metric();
  fprintf('scan %d, damage = %g, quality=%g\n', k, damage, quality);
  
end




%%
%
slice = 25:511-25;

imm_idx = 3;
mu = 100;
im_master = pixmats{imm_idx};
% im_master = SplitBregmanROF(im_master, mu, 0.001);
im_master = im_master - mean(im_master(:));
% figure, imagesc(im_master), colormap('gray')


figure(3000);clf;

[ha, pos] = tight_subplot(2, 2, [.03, .005 ], .027, .005);
thresh = (20/7)*(1/1000)*20;

fprintf('---------------------------------------------------\n');
ax_iter = 1;
for k=[1,2,3,4,7] %length(rast_exps)
  if k== imm_idx
    continue;
  end
  imk = pixmats{k};
    imk = SplitBregmanROF(imk, mu, 0.001);
  imk = imk - mean(imk(:));
  imk_slice = imk(slice, slice);
  im1_ontok_fit = norm_align(imk_slice, im_master);
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice);
  dmg = rast_exps{k}.damage_metric();
  mx = max(imk_slice(:));
  mn = min(imk_slice(:));
  stit = sprintf('PSNR=%.2f, SSIM=%.2f, Damage=%.4g', psn_1k, ssm_1k, dmg);
  
  fprintf("%s\n", stit);
  ime=imk_slice - im1_ontok_fit; %, 'parent', h(k), 'Scaling', 'joint')
  ime = ime-mean(ime(:));
  imagesc(ha(ax_iter), ime, [-2*thresh, thresh*2]); %, [-0.5*thresh, 0.5*thresh]);
  colormap('gray')
  title(ha(ax_iter), stit)
  set(ha(ax_iter), 'YTick', [], 'XTick', [])
  ax_iter = ax_iter + 1;

end



%%


function plot_raster_data(pixmat2, figbase, stit)

  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
  F10 = figure(figbase+2); clf
  ax1 = subplot(3,1,[1,2])
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
  
  figure(figbase+4)
  mesh(pixmat2)
  xlabel('x')
  ylabel('y')
  title(stit)
end
