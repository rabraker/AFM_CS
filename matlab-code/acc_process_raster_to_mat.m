clc
clear
% close all
addpath('functions')

dat_root = getdataroot();
% dat_root = '/media/labserver/acc2018-data/data/';


% -------- Constants --------
Ts = 40e-6;
% width = 5;
% nperiods = 256;
% npix = nperiods;


% dat_name_s{1} = 'raster_scan_5mic_2.50e-01Hz_out_9-17-2017-02.csv';
% dat_name_s{2} = 'raster_scan_5mic_1Hz_out_9-17-2017-02.csv';
% dat_name_s{3} = 'raster_scan_5mic_2.50e+00Hz_out_9-17-2017-01.csv';
% dat_name_s{4} = 'raster_scan_5mic_5Hz_out_9-17-2017-01.csv';
% dat_name_s{5} = 'raster_scan_5mic_10Hz_out_9-17-2017-01.csv';

% 9-17-2017
if 0
dat_name_s{1} = 'raster_scan_512pix_5mic_2.50e-01Hz_out_9-17-2017-01.csv';
dat_name_s{2} = 'raster_scan_512pix_5mic_01Hz_out_9-17-2017-01.csv';
dat_name_s{3} = 'raster_scan_512pix_5mic_2.50e+00Hz_out_9-17-2017-01.csv';
dat_name_s{4} = 'raster_scan_512pix_5mic_05Hz_out_9-17-2017-01.csv';
dat_name_s{5} = 'raster_scan_512pix_5mic_10Hz_out_9-17-2017-01.csv';
end
dat_name_s{1} = 'raster_scan_512pix_5mic_01Hz_out_9-17-2017-01.csv';
% 9-18-2017
% dat_name_s{1} = 'raster_scan_512pix_5mic_1.00e-01Hz_out_9-18-2017-03.csv';
% dat_name_s{2} = 'raster_scan_512pix_5mic_2.50e-01Hz_out_9-18-2017-01.csv';
% dat_name_s{3} = 'raster_scan_512pix_5mic_01Hz_out_9-18-2017-01.csv';
% dat_name_s{4} = 'raster_scan_512pix_5mic_2.50e+00Hz_out_9-18-2017-01.csv';
% dat_name_s{5} = 'raster_scan_512pix_5mic_05Hz_out_9-18-2017-01.csv';
% dat_name_s{6} = 'raster_scan_512pix_5mic_10Hz_out_9-18-2017-01.csv';

sub_fold = '5microns/9-22-2017';

% % 9-18-2017 20-microns
dat_name_s = {};
% dat_name_s{1} = 'raster_scan_512pix_20mic_1.00e-01Hz_out_9-23-2017-03.csv';
dat_name_s{1} = 'raster_scan_512pix_20mic_05Hz_out_9-23-2017-01.csv';
% dat_name_s{3} = 'raster_scan_512pix_20mic_2.50e+00Hz_out_9-23-2017-01.csv';
% % dat_name_s{5} = 'raster_scan_512pix_20mic_1.00e-01Hz_out_9-18-2017-03.csv';
sub_fold = '20microns/9-23-2017';

% dat_name_s = {};
% dat_name_s{1} = 'raster_scan_512pix_5mic_2.50e-01Hz_out_9-23-2017-02.csv';
% dat_name_s{2} = 'raster_scan_512pix_5mic_01Hz_out_9-23-2017-04.csv';
% dat_name_s{3} = 'raster_scan_512pix_5mic_05Hz_out_9-23-2017-03.csv';
% dat_name_s{4} = 'raster_scan_512pix_5mic_10Hz_out_9-23-2017-01.csv';




freqs = zeros(1,length(dat_name_s));

verbose = 2;
raster_dat_s = [];
for i=1:length(dat_name_s)
%     dat_path= fullfile(dat_root,'raster/9-14-2017', dat_name_s{i});
    dat_path= fullfile(dat_root,'raster', sub_fold, dat_name_s{i});
    parent_name = get_parent_name(dat_name_s{i}, '_out_');
    parent_path = fullfile(dat_root, 'raster', sub_fold, parent_name);
    parent_meta_path = strrep(parent_path, '.csv', '.mat');
    meta_path = strrep(dat_path, '.csv', '-meta.mat');

    fprintf('File: %s\n', dat_name_s{i}); 
%     load(meta_path)
%     keyboard
    mat_path = strrep(dat_path, '.csv', '.mat');
    if exist(mat_path, 'file') ~=2 
%         keyboard
        fprintf('Reading in raw data file...\n');
        rdat.datmat = csvread(dat_path);
        rdat.parent_dat = csvread(parent_path);
        load(parent_meta_path);
        rdat.meta = meta;
        freqs(i) = rdat.meta.raster_freq;

        xyref = reshape(rdat.parent_dat, 2, [])';
        rdat.xref = xyref(:,1);
        rdat.Ts = Ts;
        rdat.npix = rdat.meta.npix;
        volts2micron = 50/10;
        micron2pix = rdat.npix/rdat.meta.width;
        rdat.volts2pix = volts2micron * micron2pix;
        rdat.width = rdat.meta.width;

        rdat.samps_per_period = size(rdat.parent_dat,1)/2; % twice as many in here for x & y.
        rdat.samps_per_line = rdat.samps_per_period/2;
        fprintf('Processing raster data...\n')
        rdat.datmat = rdat.datmat([1:rdat.npix*rdat.samps_per_period], :);

        [pixmat, pixelifsampled, m_s] = bin_raster_really_slow(rdat.datmat(:,[1,2,4]),...
                                    rdat.npix, rdat.samps_per_period, rdat.volts2pix);
        I_fit = detrend_plane(pixmat);
%         I_fit = detrend_sampled_plane(pixmat, pixelifsampled);
%         I_fit = (I_fit - min(min(I_fit))).*pixelifsampled;
        rdat.pixmat = pixmat;
        rdat.I_fit = I_fit;
        rdat.pixelifsampled = pixelifsampled;
        if verbose
            figure;
            imshow(I_fit, [min(min(I_fit)), max(max(I_fit))])
            drawnow
        end

        rdat.freq = freqs(i);
        raster_dat_s = [raster_dat_s, rdat];
        fprintf('saving raster data...\n');
        rdat.datmat=[];
         save(mat_path, 'rdat')
    else
        load(mat_path);
        fprintf('File already processed. Skipping...\n')
        raster_dat_s = [raster_dat_s, rdat];
        freqs(i) = rdat.meta.raster_freq;
    end

end



%%
% close all
master = raster_dat_s(1);
max_s = zeros(1,4);
width = 5;
height = width;
xdata = [0, width];
ydata = [0, height];

clc
thresh = (20/7)*(1/1000)*120*2;
for i=1:length(freqs)
    
    I_fit = raster_dat_s(i).I_fit;
    [kmin, kmax] = find_raster_extents(raster_dat_s(i));
    I_temp = I_fit(:,kmin:kmax);
    I_fit = I_fit-mean(I_temp(:));
    I_fit = max(I_fit, -thresh);
    I_fit = min(I_fit, thresh);
%     I_fit = 
    I_temp = I_temp-mean(I_temp(:));
    I_temp = I_fit(:,kmin:kmax);
    [min_s(i), ind_min] = min(min(I_temp));
     max_s(i) = max(max(I_fit));
   %    I_fit = I_fit - min_s(i);
%    I_temp = I_temp - min_s(i);
   % set the unscanned area to the level of the substrate:
%    I_fit(:,1:kmin-1) = max(max(I_temp));
%    I_fit(:,kmax+1:end) = max(max(I_temp));
    
    raster_dat_s(i).I_fit2 = I_fit;

end

maxx = max(max_s);
minn = min(min_s);


% cs-traj-10perc-500nm-5mic-01Hz_out_9-15-2017-02_img-data.mat
cs_root =  fullfile(dat_root,'cs-data');
cs_name_s{1} = 'cs-traj-512pix-8perc-500nm-5mic-5.00e-01Hz_out_9-17-2017-01.csv';
cs_name_s{2} = 'cs-traj-512pix-10perc-500nm-5mic-5.00e-01Hz_out_9-17-2017-01.csv';

clc
% cs_name_s = {};
% % cs_name_s{1} = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_9-18-2017-02.csv';
% % cs_name_s{2} = 'cs-traj-512pix-15perc-500nm-5mic-01Hz_out_9-18-2017-01.csv';

% cs_name_s{3} = 'cs-traj-512pix-10perc-500nm-5mic-5.00e-01Hz_out_9-18-2017-02.csv'
% ca_name_s{3} = 'cs-traj-512pix-15perc-500nm-5mic-01Hz_out_9-18-2017-01.csv'

% cs_name_s{1} = 'cs-traj-512pix-10perc-2000nm-20mic-01Hz_out_9-18-2017-01.csv';
% cs_name_s{2} = 'cs-traj-512pix-3perc-2000nm-20mic-01Hz_out_9-18-2017-01.csv';

max_cs = zeros(1,length(cs_name_s));
min_cs = zeros(1,length(cs_name_s));
for i=1:length(cs_name_s)
    img_data_name = strrep(cs_name_s{i}, '.csv', '_img-data.mat');
    parent_name = get_parent_name(cs_name_s{i}, '_out_');
    meta_name = strrep(parent_name, '.csv', '.mat');
    load(fullfile(cs_root,sub_fold, meta_name));
    
    load(fullfile(cs_root, sub_fold, img_data_name))
    img_data.CsExpMetaIn = CsExpMetaIn;
    pixelifsampled = img_data.CsExpMetaIn.pixelifsampled;
    pix = img_data.CsExpMetaIn.npix;
    actual_sub_sample_frac = length(find(pixelifsampled == 1))/pix^2;
    
    cs_img_data_s{i} = img_data;
    cs_img_data_s{i}.img_data.perc = actual_sub_sample_frac*100;
    bp_im = img_data.bp_im;
    
    bp_im = bp_im - mean(bp_im(:));
%     bp_im = bp_im - min(bp_im(:));
    min_cs(i) = min(min(bp_im));
    max_cs(i) = max(max(bp_im));
    
    cs_img_data_s{i}.bp_im = bp_im;
end

DRng = max(max(max_s), max(max_cs)) - min(min(min_cs), min(min_s));

for i=1:length(freqs)
    
    raster_dat_s(i).I_fit_normalized = raster_dat_s(i).I_fit2;
    
    figure(i);
    %   imshow(I_fit, [min(min(I_fit)), max(max(I_fit))], 'XData', xdata, 'YData', ydata)
    imshow(raster_dat_s(i).I_fit_normalized, [min_s(i), max_s(i)])
%     , 'XData', xdata, 'YData', ydata)
    axis('on') 
    drawnow
    stit = sprintf('%2d Hz', freqs(i));
    title(stit)
    ax = gca;

    grid on
    set(ax,'xminorgrid','on','yminorgrid','on')
    colorbar

    
end


for i=1:length(cs_name_s)
    cs_img_data_s{i}.bp_im_normalized = cs_img_data_s{i}.bp_im ;
    figure(i+length(freqs));
    imshow(cs_img_data_s{i}.bp_im_normalized, [min_cs(i), max_cs(i)])
%     [0, DRng])
%     ,  'XData', xdata, 'YData', ydata)
    axis('on') 
    drawnow
    title(sprintf('bp: %d%%', floor(cs_img_data_s{i}.img_data.perc)))
    ax = gca;

    grid on
      set(ax,'xminorgrid','on','yminorgrid','on')
      colorbar
end
%%
clc
ssim_s = zeros(1, length(freqs)-1);
psnr_s = zeros(1, length(freqs)-1);

st_head = sprintf('-----|');
st_ssm  = sprintf('ssim |');
st_psnr = sprintf('psnr |');
k = 1;
L = 1;
REF = raster_dat_s(1).I_fit_normalized(1:end-L+1,1:end-k+1);
expns = [1,1,1];
figure(1000);
subplot(2,1,2)
hold on
for iter=1:length(freqs)
    
%     X = cs_img_data_s{2}.bp_im;
    
    Y = raster_dat_s(iter).I_fit_normalized(L:end,k:end);

%     mn = min(min(REF(:)), min(Y(:)));
%     mx = max(max(REF(:)), max(Y(:)));
     mn = min(REF(:));
     mx = max(REF(:));
    DRng_i = mx-mn;

%    ssim_s(iter) = mean(ssms);
%    if ssim_s(iter) <0
%        keyboard
%    end
   ssim_s(iter) = ssim(Y, REF, 'DynamicRange', DRng);
   psnr_s(iter) = psnr(Y, REF, DRng_i);
   if isinf(psnr_s(iter))
       st_fill = '   ';
   else
       st_fill = '';
   end
   st_head = sprintf('%s  %.3f Hz|', st_head, freqs(iter));
   st_ssm = sprintf('%s  %.3f   |', st_ssm, ssim_s(iter));
   st_psnr = sprintf('%s  %.3f  %s|', st_psnr, psnr_s(iter),st_fill);
   plot(Y(40,:))
end
%
% legend('0.1Hz', '0.25Hz', '1Hz', '2.5Hz')
legend( '0.25Hz', '1Hz', '2.5Hz', '5Hz', '10Hz')
% k = 23
k = 25;
L = 1;
REF = raster_dat_s(1).I_fit_normalized(1:end-L+1,k:end);
for i=1:length(cs_img_data_s)
    Y = cs_img_data_s{i}.bp_im_normalized(L:end, 1:end-k+1);
%     mn = min(min(REF(:)), min(Y(:)));
%     mx = max(max(REF(:)), max(Y(:)));
     mn = min(REF(:));
     mx = max(REF(:));
    DRng_i = mx-mn;
    
    st_head = sprintf('%s bp:%.2f%%   |', st_head, cs_img_data_s{i}.img_data.perc);
   
    ssim_bp = ssim(Y, REF, 'DynamicRange', DRng_i);
%     ssim_bp = mean(ssms);
    psnr_bp = psnr(Y, REF, DRng_i);
    st_ssm = sprintf('%s    %.3g |', st_ssm, ssim_bp);

    st_psnr = sprintf('%s    %.3g |', st_psnr, psnr_bp);

end

st_result = sprintf('%s\n%s\n%s\n', st_head, st_ssm, st_psnr)

%%
% ------ Use my muPathMask
REF = raster_dat_s(1).I_fit_normalized;
REF_samp = REF.*cs_img_data_s{1}.pixelifsampled;

bp_im_mine = bp_reconstruct(REF_samp, cs_img_data_s{1}.pixelifsampled)

% ---- Use Yufans mu-path mask
pixelifsampled_yufan = muPathMaskGen(53, 512, 512, 0.083263397216797, false);
REF_samp_yufan = REF.*pixelifsampled_yufan;


bp_im_yufan = bp_reconstruct(REF_samp_yufan, pixelifsampled_yufan);
%%
figure(1)
subplot(2,3,1)
imagesc(REF)
colormap gray
title('ground truth image')


subplot(2,3,2)
imagesc(cs_img_data_s{1}.pixelifsampled); 
colormap gray;
title('My mu-path mask')
subplot(2,3,3)
imagesc(bp_im_mine)
colormap gray
title('BP reconstruction, my-mask')

subplot(2,3,4)
imagesc(cs_img_data_s{1}.bp_im_normalized)
colormap gray
title('BP Reconstruction, AFM')




subplot(2,3,5)
imagesc(pixelifsampled_yufan)
colormap gray
title('Yufans mu-path mask')

subplot(2,3,6)
imagesc(bp_im_yufan)
title('BP reconstruction, Yufans mask')





%%
% close all
im1 = raster_dat_s(2).I_fit_normalized;
im2 = raster_dat_s(3).I_fit_normalized;


x=15; y=15;
REF=256-x; Y=256-y;
im2_sub = im2(y:Y,x:REF);

figure(1)
subplot(2,3,1)
imagesc(flipud(im1))
colormap gray
title('master')
hold on
plot([x,x,REF,REF,x], [y, Y,Y,y,y], 'r')

subplot(2,3,2)
imagesc(flipud(im2))
colormap gray
title('drifted')
hold on
plot([x,x,REF,REF,x], [y, Y,Y,y,y], 'r')

subplot(2,3,3)

imagesc(flipud(im2_sub))
colormap gray


nim = im1 - mean(mean(im1));
n_im2 = im2_sub - mean(mean(im2_sub));

ccr = xcorr2(nim, n_im2);

[ssr, snd] = max(ccr(:));
[ij, ji] = ind2sub(size(ccr), snd)
subplot(2,3,4)
plot(ccr(:))
title('cross correlation')
hold on
plot(snd,ssr,'or')
hold off
text(snd*1.05,ssr,'Maximum')



im1_new = im1(x:REF, y:Y);

% shft=0
im2_new = ones(256,256);
im2_new(ij:-1:ij-size(im2_sub,1)+1,ji:-1:ji-size(im2_sub,2)+1) = rot90(im2_sub,2);

% im2_new = im1(ij-size(im2_sub,1):ij,ji-size(im2_sub,2):ij);

subplot(2,3,5)
imagesc(flipud(im1_new))
ax=gca;
% axis image off
colormap gray
title('master sliced')
axis('on')
grid on
set(ax,'xminorgrid','on','yminorgrid','on')

subplot(2,3,6)
imagesc(flipud(im2_new))
% axis image off
colormap gray
title('shifted slice')
ax = gca;
axis('on')
grid on
set(ax,'xminorgrid','on','yminorgrid','on')

