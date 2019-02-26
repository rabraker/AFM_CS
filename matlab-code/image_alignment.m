% This script reads in the full, raw data from a raster scan run. This data
% will be produced by the vi play-raster-scan.vi. That vi produces two .csv
% files. One contains the the pre-processed data that labview does in
% realtime for visualization purposed. For slow scans, using that data is
% sufficient. The vi also produces a csv file with all of the raw data. For
% faster scans, we want to use that data, so we can process it better.
% That's what this script does.
clc


init_paths();
import scanning_v1.*

clear size
size_dir = '5microns';

% dat_root = PATHS.raster_image_data(size_dir, '2-25-2019');
% dat_name1 = 'raster_scan_256pix_5mic_01Hz_out_2-25-2019-01.csv';
% dat_name2 = 'raster_scan_256pix_5mic_01Hz_out_2-25-2019-02.csv';
dat_root = PATHS.raster_image_data(size_dir, '2-25-2019');


dat_name1 = 'raster_scan_256pix_5mic_01Hz_out_2-25-2019-03.csv';
dat_name2 = 'raster_scan_256pix_5mic_01Hz_out_2-25-2019-04.csv';
dat_names = {dat_name1, dat_name2};

pixmaps = cell(length(dat_names), 1);
for k=1:length(dat_names)
  
  raster_paths = get_raster_paths(dat_root, dat_names{k});

  rast_exps{k} = RasterExp(raster_paths);
end
%%
npix = rast_exps{1}.npix;
Ts = rast_exps{1}.Ts;
width = rast_exps{1}.width;
stit = sprintf('Scan %d', k)

x1s = [5, 18];
x2s = [235 253];
figbase = 10;
for k=1:2
  rast_exp2 = copy(rast_exps{k});
  wo = 2*pi*rast_exp2.meta_data.raster_freq;
%   uz = detrend(rast_exp2.uz);
  
%   y_est = fourier_tri(uz, wo);
%   rast_exp2.uz = real(fft_highpass(uz, Ts, .05));
%   rast_exp2.uz = rast_exp2.uz- y_est;
  % rast_exp.bin_raster_really_slow();
  
  rast_exp2.bin_raster_really_slow();
  
  pixmats{k} = rast_exp2.pix_mat(1:end-1, 1:end-1);
  pixmats{k} = pin_along_column(pixmats{k}, x1s(k), x2s(k));
  pixmats{k} = pixmats{k} - mean(pixmats{k}(:));
  stit = sprintf('scan %d', k);
  plot_data(pixmats{k}, figbase*k, stit)
  
  
end

%%


im1 = pixmats{1};
im2 = pixmats{2};
im2 = im2(10:256-10, 10:256-10);

im1_fit = norm_align(im2, im1);

figure(11); clf
imshowpair(im1_fit, im2)    
title('fit')
%%
figure(12)
imshowpair(pixmats{1}, pixmats{2})
title('original')
%%
figure(6)
imagesc(im1_fit, [-0.05, 0.05])
colormap('gray')
figure(7)
imagesc(im2, [-0.05, 0.05])
colormap('gray')

[psn, ssm] = ssim_psnr_norm(pixmats{1}(10:256-10, 10:256-10), pixmats{2}(10:256-10, 10:256-10))
[psn, ssm] = ssim_psnr_norm(im1_fit, im2);

fprintf('psnr: %.4f\n', psn);
fprintf('ssm: %.4f\n', ssm);
%%

rast_exp2 = rast_exp;
rast_exp2.uz = detrend(rast_exp2.uz);

wo = 2*pi*rast_exp2.meta_data.raster_freq

% y_est = fourier_tri(rast_exp2.uz, wo);
% rast_exp2.uz = rast_exp2.uz - y_est;

uu_filt = real(fft_notch(rast_exp2.uz, AFM.Ts, 0.49, 0.51));
uu_filt = real(fft_notch(uu_filt, AFM.Ts, 0.01, 0.1));
rast_exp2.uz = uu_filt;

plot_data(rast_exp2, figbase*k, stit);


function pixmat = pin_along_column(pixmat, x1, x2)

  xs = 1:size(pixmat,2);
  for k=1:size(pixmat,1)
    y1 = pixmat(k, x1);
    y2 = pixmat(k, x2);
    
    m = (y2-y1)/(x2-x1);
    b = y1 - m*x1;
    line = m*xs + b;
    pixmat(k,:) = pixmat(k,:) - line;
  end
  

end

function plot_data(pixmat2, figbase, stit)

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
  grid(ax1, 'on')
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

function y_est = fourier_tri(uu, wo)
  N = length(uu);
  t = double(0:1:N-1)'*(AFM.Ts);
%   plot(t, uu)
  
  
  num_coef = 1;
  clc
    
  ao_s = []
  bo_s = []
  y_est = 0;
  for k=1:num_coef
    
    ys = sin(k*wo*t);
    yc = cos(k*wo*t);
    
    ao = ys'*uu*(2/N)
    bo = yc'*uu*(2/N)
    
    ao_s(k) = ao;
    bo_s(k) = bo;
    
    y_est = y_est + ao * ys + bo*yc;
  end
end
%%
% % % optimizer = registration.optimizer.RegularStepGradientDescent;
% % % metric = registration.metric.MeanSquares;
% % [optimizer, metric] = imregconfig('Multimodal');
% % [pixmat_reg, reg_dat] = imregister(pixmats{1}, pixmats{2}, 'rigid', optimizer, metric);
% % 
% % % imshowpair(pixmats{1}, pixmats{2})
% % % pixmat_reg = imregister(pixmats{1}, pixmats{2}, 'Similarity', optimizer, metric);
% % 
% % figure(6)
% % imagesc(pixmat_reg-mean(pixmat_reg(:)), [-0.05, 0.05])
% % colormap('gray')
% % figure(7)
% % imagesc(pixmats{2}-mean(pixmats{2}(:)), [-0.05, 0.05]);
% % colormap('gray')
% % 
% % figure(8)
% % imshowpair(pixmat_reg, pixmats{2})

