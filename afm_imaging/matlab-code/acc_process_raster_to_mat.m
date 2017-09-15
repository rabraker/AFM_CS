clc
clear
close all
addpath('functions')

% dat_root = getdataroot();
dat_root = '/media/labserver/acc2018-data/data/';


% -------- Constants --------
Ts = 40e-6;
width = 5;
nperiods = 256;
npix = nperiods;
volts2micron = 50/10;
micron2pix = npix/width;
volts2pix = volts2micron * micron2pix;


% ------------------- 1Hz 
% % dat_name_s{1} = 'raster_scan_5mic_1Hz_out_8-26-2017-01-full.csv';
% dat_name_s{1} = 'raster_scan_5mic_1Hz_out_9-14-2017-01-full.csv';
% % ------------------- 5Hz 
% dat_name_s{2} = 'raster_scan_5mic_5Hz_out_9-14-2017-01-full.csv';
% % ------------------- 10Hz 
% dat_name_s{3} = 'raster_scan_5mic_10Hz_out_9-14-2017-01-full.csv';

if 1
    % ------------------- 1Hz 
    % dat_name_s{1} = 'raster_scan_5mic_1Hz_out_8-26-2017-01-full.csv';
    % parent_name_s{1} = 'raster_scan_5mic_1Hz.csv';
    dat_name_s{1} =  'raster_scan_5mic_5.00e-01Hz_out_9-14-2017-03.csv';
    % ------------------- 1Hz 
    dat_name_s{2} = 'raster_scan_5mic_1Hz_out_9-14-2017-04.csv';    
    % ------------------- 5Hz 
    dat_name_s{3} = 'raster_scan_5mic_5Hz_out_9-14-2017-02.csv';
    % ------------------- 10Hz 
    dat_name_s{4} = 'raster_scan_5mic_10Hz_out_9-14-2017-02.csv';
end



freqs = [.5, 1, 5, 10];
verbose = 2;
raster_dat_s = [];
for i=1:length(freqs)
%     dat_path= fullfile(dat_root,'raster/9-14-2017', dat_name_s{i});
    dat_path= fullfile(dat_root,'raster/9-14-2017_2', dat_name_s{i});
    parent_name = get_parent_name(dat_name_s{i}, '_out_');
    parent_path = fullfile(dat_root, 'raster/9-14-2017_2', parent_name);
    
    mat_path = strrep(dat_path, '.csv', '.mat');
        if exist(mat_path, 'file') ~=2

        rdat.datmat = csvread(dat_path);
%         keyboard
        rdat.parent_dat = csvread(parent_path);

        xyref = reshape(rdat.parent_dat, 2, [])';
        rdat.xref = xyref(:,1);
        rdat.Ts = Ts;
        rdat.npix = npix;
        rdat.volts2pix = volts2pix;
        rdat.width = width;

        rdat.samps_per_period = size(rdat.parent_dat,1)/2; % twice as many in here for x & y.
        rdat.samps_per_line = rdat.samps_per_period/2;
        rdat.datmat = rdat.datmat([1:rdat.npix*rdat.samps_per_period], :);

        [pixmat, pixelifsampled] = bin_raster_slow(rdat.datmat(:,[1,2,4]),...
                                    rdat.npix, rdat.samps_per_period, rdat.volts2pix);
        I_fit = detrend_sampled_plane(pixmat, pixelifsampled);
        I_fit = (I_fit - min(min(I_fit))).*pixelifsampled;
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

        save(mat_path, 'rdat')
    else
        load(mat_path);
        fprintf('File: %s\nalready processed. Skipping...\n', dat_name_s{i}); 
        raster_dat_s = [raster_dat_s, rdat];
    end

end
%%
close all
master = raster_dat_s(1);
max_s = zeros(1,4);
for i=1:4
    I_fit = raster_dat_s(i).I_fit;
    [kmin, kmax] = find_raster_extents(raster_dat_s(i));
    I_temp = I_fit(:,kmin:kmax);
   min_s(i) = min(min(I_temp))
   I_fit = I_fit - min_s(i);
   
   max_s(i) = max(max(I_fit));
   
  figure(i);
  imshow(I_fit, [min(min(I_fit)), max(max(I_fit))])
  drawnow
end
%%
maxx = max(max_s);

for i=1:4
    
    I_fit = raster_dat_s(i).I_fit;
    %scale to 255;
    I_fit = I_fit*(255/maxx);
    raster_dat_s(i).I_fit = I_fit;
    figure(i)
    imshow(I_fit, [0, 255])
end
%%
clc
for i=2:4
    X = raster_dat_s(1).I_fit;
    Y = raster_dat_s(i).I_fit;
    
   ssim(X, Y)
   psnr_xy(X, Y)
end















