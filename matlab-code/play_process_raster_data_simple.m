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

dat_root = PATHS.raster_image_data(size_dir, '2-19-2019');
dat_name1 = 'raster_scan_128pix_5mic_01Hz_out_2-19-2019-03.csv';
dat_name2 = 'raster_scan_128pix_5mic_01Hz_out_2-19-2019-04.csv';
dat_names = {dat_name1, dat_name2};

pixmaps = cell(2)
for k=1:2
  
raster_paths = get_raster_paths(dat_root, dat_names{k});

rast_exp = RasterExp(raster_paths);

npix = rast_exp.npix;
Ts = rast_exp.Ts;
width = rast_exp.width;
stit = sprintf('Scan %d', k)


pixmaps{k} = plot_data(rast_exp, figbase*k, stit)
end
%%

pixmap1 = pixmaps{1};
pixmap1 = pixmap1 - min(pixmap1(:));
pixmap2 = pixmaps{2};
pixmap2 = pixmap2 - min(pixmap2(:));

mx1 = max(abs(pixmap1(:)));
mx2 = max(abs(pixmap2(:)));
mx = max(mx1, mx2);

psn = psnr(pixmap1, pixmap2, mx);
ssm = ssim(pixmap1, pixmap2, 'DynamicRange', mx);

fprintf('psnr: %.4f\n', ssm);
fprintf('ssm: %.4f\n', psn);
%%


function pixmat2 = plot_data(rast_exp, figbase, stit)
clc
% rast_exp.bin_raster_really_slow(@detrend);
rast_exp.bin_raster_really_slow();

pixmat2 = rast_exp.pix_mat;

thresh = (20/7)*(1/1000)*20;
pixmat2 = pixmat2 - mean(pixmat2(:));
F10 = figure(figbase+2); clf
ax1 = gca();

f11 = figure(figbase+3); clf
ax2 = gca();

lo = min(min(pixmat2));
hi = max(max(pixmat2));


pixmat2 = detrend2(detrend2(pixmat2(10:end-10,10:end-25)));
imshow_dataview(flipud(pixmat2 - mean(pixmat2(:))), [-thresh, thresh], ax1, ax2)
grid(ax1, 'on')
colormap(ax1, 'parula')
grid(ax2, 'on')
ax1.GridAlpha = 1;
ax2.GridAlpha = 1;
title(ax1, stit)
title(ax2, stit)

figure(12+figbase)
mesh(pixmat2)
xlabel('x')
ylabel('y')
title(stit)
end
%%
% x = [1:512];
% y = x;
% [xx, yy, zz] = prepareSurfaceData(x, y, pixmat2);
% %%
% pixmat3 = pixmat2*0;
% f= fit([xx, yy], zz, 'poly23')
% for k=1:512
%   pixmat3(k, :) = pixmat2(k,:) - f(k, x');
% end
% [xx, yy, zzz] = prepareSurfaceData(x, y, pixmat3);
% figure, plot([xx, yy], zzz)


