clc
clear
addpath functions/scanning_v1/
% initialize paths.
init_paths();

Ts = 40e-6;

size_dir = '5microns';

data_root = PATHS.cs_image_data(size_dir, '3-3-2019');

cs_exp_data_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-isconnect_out_3-3-2019-01.csv';
cs_exp_data_name_s{2} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz-250prescan-notconnect_out_3-3-2019-01.csv';


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

cs_paths = get_cs_paths(data_root, cs_exp_data_name_s{2});

hole_depth = (20);

cs_exp = CsExp(cs_paths, 'feature_height', hole_depth);% 'gg', gg);
cs_exp.uz = detrend(cs_exp.uz);
cs_exp.print_state_times();

figbase = 20;

%
[~, axs] = make_traj_figs(figbase);
cs_exp.plot_all_cycles(axs{:});



cs_exp.process_cs_data(false, []);
pixelifsampled = cs_exp.pix_mask;
I = cs_exp.Img_raw;
fprintf('finished processing raw CS data...\n');

fprintf('nperc=%.3f\n', sum(cs_exp.pix_mask(:))/cs_exp.npix^2);

ht = cs_exp.feature_height;
figure(10+figbase)
ax = gca();
figure(11+figbase)
axx = gca();
imshow_dataview(cs_exp.Img_raw-mean(cs_exp.Img_raw(:)), [-ht, ht], ax, axx)

%
bp = true;

clear CsExp
clear l1qc_dct_mex
 
ht = cs_exp.feature_height;

addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib
addpath cs_simulations/splitBregmanROF_mex/
figure(13)
ax4 = gca();

if bp
  %%
  
  clear CsTools
    cs_exp.solve_bp(true, true);
  %%
  f10=figure(10); clf
  ax4 = subplot(3,1,[1,2])
  ax4_2 =subplot(3,1,3);

  ImshowDataView.setup(f10);
  cb_exp =  @(event_obj)cs_exp.dataview_callback(event_obj, ax4, ax4_2);
  ImshowDataView.imshow(cs_exp.Img_bp, [], ax4, ax4_2, cb_exp)
  title(ax4, 'BP reconstruction');
end
%%
f10=figure(10); clf
ax4 = subplot(3,1,[1,2])
ax4_2 =subplot(3,1,3);
% imshow_dataview(cs_exp.Img_bp, [-ht, ht], ax4, ax4_2)
%
im = cs_exp.Img_bp;
for k=1:size(im,1)
  im(k, :) = detrend(im(k,:));
end
ImshowDataView.setup(f10);
cb_exp =  @(event_obj)cs_exp.dataview_callback(event_obj, ax4, ax4_2);
ImshowDataView.imshow(cs_exp.Img_bp, [], ax4, ax4_2, cb_exp)
%%
dat_root = PATHS.raster_image_data(size_dir, '2-25-2019');


dat_name1 = 'raster_scan_256pix_5mic_01Hz_out_2-25-2019-03.csv';
dat_name2 = 'raster_scan_256pix_5mic_01Hz_out_2-25-2019-04.csv';
dat_names = {dat_name1, dat_name2};

pixmaps = cell(length(dat_names), 1);
for k=1:length(dat_names)
  
  raster_paths = get_raster_paths(dat_root, dat_names{k});

  rast_exps{k} = RasterExp(raster_paths);
end

npix = rast_exps{1}.npix;
Ts = rast_exps{1}.Ts;
width = rast_exps{1}.width;
stit = sprintf('Scan %d', k)
%%
clc
x1s = [23, 24];
x2s = [235 253];
figbase = 10;
for k=1:2
  rast_exp2 = copy(rast_exps{k});
  wo = 2*pi*rast_exp2.meta_data.raster_freq;
%   uz = detrend(rast_exp2.uz);

  rast_exp2.bin_raster_really_slow();
  
  pixmats_raw{k} = rast_exp2.pix_mat(1:end-1, 1:end-1);
  pixmats{k} = pin_along_column(pixmats_raw{k}, x1s(k), x2s(k));
  pixmats{k} = pixmats{k} - mean(pixmats{k}(:));
  stit = sprintf('(pinned) scan %d', k);
  plot_raster_data(pixmats{k}, figbase*k, stit)
%   stit = sprintf('(raw) scan %d', k);
%   plot_raster_data(pixmats_raw{k}, (figbase-5)*k, stit)
  
end
%%
Img_filt = cs_exp.Img_bp - min(cs_exp.Img_bp(:));
% Img_filt = im - min(im(:)); %min(cs_exp.Img_bp(:));
% Img_filt = Img_filt*255/max(Img_filt(:));
% for k=1:256
%   Img_filt(k,:) = detrend(Img_filt(k,:));
% end


Img_filt = SplitBregmanROF(Img_filt, 30, 0.001);
figure(11)
% Img_filt = pin_along_column(Img_filt, 49, 464);
mesh(Img_filt), colormap('gray')

%%

slice = 10:256-10;
im1 = pixmats{1} - mean(pixmats{1}(:));
im1 = SplitBregmanROF(im1, 45, 0.001);
im2 = pixmats{2} - mean(pixmats{2}(:));
im2 = SplitBregmanROF(im2, 45, 0.001);
im2_slice = im2(slice, slice);

im3 = Img_filt(slice, slice);
im3 = im3 - mean(im3(:));

im1_onto2_fit = norm_align(im2_slice, im1);
im1_onto3_fit = norm_align(im3, im1);

figure(11); clf
imshowpair(im1_onto2_fit, im2_slice)   
title('fit (raster)')
[psn_12, ssm_12] = ssim_psnr_norm(im1_onto2_fit, im2_slice);


figure(12); clf
imshowpair(im1_onto3_fit, im3); 
title('fit (CS)')
[psn_13, ssm_13] = ssim_psnr_norm(im1_onto3_fit, im3);

fprintf('---------------------------------------------------\n');
fprintf('(raster) psnr: %.4f, ssm: %.4f\n  ', psn_12, ssm_12);
fprintf('(CS) psnr: %.4f, ssm: %.4f\n  ', psn_13, ssm_13);


figure(51)
mesh(im1)
colormap('gray')
title('raster-1, TV')
xlabel('x')
ylabel('y')

figure(52)
mesh(im2)
colormap('gray')
title('raster-2, TV')
xlabel('x')
ylabel('y')

figure(53)
mesh(im3)
colormap('gray')
title('CS, TV')
xlabel('x')
ylabel('y')
%%
figure(61)
mesh(pixmats{1})
% colormap('gray')
title('raster-1, raw')

figure(62)
mesh(pixmats{2})
% colormap('gray')
title('raster-2, raw')

figure(63)
mesh(cs_exp.Img_bp)
% colormap('gray')
title('CS, TV')
%%
function plot_raster_data(pixmat2, figbase, stit)

  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
  F10 = figure(figbase+2); clf
  ax1 = gca();
  
  f11 = figure(figbase+3); clf
  ax2 = gca();
  
  lo = min(min(pixmat2));
  hi = max(max(pixmat2));
  
  
  % pixmat2 = detrend2(detrend2(pixmat2(10:end-10,10:end-25)));
  
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

function [figs, axs] = make_traj_figs(figbase)
  Fig_uz = figure(20+figbase); clf
  
  ax1 = gca();
  Fig_ze = figure(30+figbase); clf
  ax2 = gca();
  
  Fig_x = figure(40+figbase); clf
  ax3 = gca();
  Fig_y = figure(50+figbase); clf
  ax4 = gca();
  
  figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
  axs = {ax1, ax2, ax3, ax4};
end