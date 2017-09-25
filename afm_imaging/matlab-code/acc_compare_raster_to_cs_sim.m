clear
clc

addpath('functions')

cs_sim_root = '/media/labserver/data/cs-sim-data';
csin_data_root = '/media/labserver/data/cs-data';
raster_root = '/media/labserver/data/raster';
% raster_root = '/home/arnold/gradschool/afm-cs/afm_imaging/acc2018-data/raster';


% csin_data_name1 = 'cs-traj-7perc-500nm-5mic-01Hz.csv';
% csin_data_name2 = 'cs-traj-10perc-500nm-5mic-01Hz.csv';
% csin_data_name3 = 'cs-traj-15perc-500nm-5mic-01Hz.csv';


% raster_name = 'raster_scan_5mic_1Hz_out_8-26-2017-01-full.csv';
% raster_parent_name = 'raster_scan_5mic_1Hz.csv';


raster_parent_name = 'raster_scan__256pix_10mic_5.00e-01Hz.mat';
raster_parent_meta = 'raster_scan__256pix_10mic_5.00e-01Hz.csv';

raster_name = 'raster_scan__256pix_10mic_5.00e-02Hz_out_9-10-2017-02-full.csv'

csin_data_name1 = 'cs-traj-7perc-500nm-5mic-01Hz.csv';
csin_data_name2 = 'cs-traj-10perc-500nm-5mic-01Hz.csv';
csin_data_name3 = 'cs-traj-15perc-500nm-5mic-01Hz.csv';


% raster_name = 'raster_scan_5mic_1Hz_out_8-26-2017-01-full.csv';
% raster_parent_name = 'raster_scan_5mic_1Hz.csv';



csin_data_path_s{1} = fullfile(csin_data_root, csin_data_name1);
csin_data_path_s{2} = fullfile(csin_data_root, csin_data_name2);
csin_data_path_s{3} = fullfile(csin_data_root, csin_data_name3);



raster_path = fullfile(raster_root, raster_name);
raster_parent_path = fullfile(raster_root, raster_parent_name);

% raster_data = csvread(raster_path);

% cs-traj15-500-1Hz.csv
% -------- Constants --------
rdat.Ts = 40e-6;
rdat.width = 10;

rdat.npix = 256;

% microns_per_volt = 50/10;
micron2pix = rdat.npix/rdat.width;
rdat.volts2pix = volts2microns * micron2pix;
rdat.pix_per_volt = (rdat.npix/rdat.width)*volts2microns;

rdat = raster_full_to_img(raster_path, raster_parent_path, rdat);

figure(1)
imshow(rdat.I_fit, [min(min(rdat.I_fit)), max(max(rdat.I_fit))]);


[kmin, kmax] = find_raster_extents(rdat);

rdat.I_fit = rdat.I_fit(:, kmin:kmax);
rdat.I_fit = detrend_plane(rdat.I_fit);

%%
figure(2)
subplot(2,2,1)
imshow_sane(rdat.I_fit, gca, rdat.width, rdat.width);
title('1 Hz raster')

%%

cs_sim_data.csin_data_paths = csin_data_path_s;
cs_sim_data.bp_im_s = {};
cs_sim_data.pixelifsampled_s = {};
cs_sim_data.rdat = rdat;

for k=1:length(csin_data_path_s)
    
    
    csin_data_path_s{k}
    
    [bp_im, pixelifsampled] = cs_sim(rdat, csin_data_path_s{k});
    
    
    
    
    perc = sum(sum(pixelifsampled))/rdat.npix^2;
    figure(2)
    subplot(2,2,k+1)
    imshow(bp_im, [min(min(bp_im)), max(max(bp_im))]);
    
    stit = sprintf('%.2d %%', perc);
    title(stit)
    
    
    figure(100+k)
    imshow_sane(pixelifsampled, gca, rdat.width, rdat.width);
    title(stit)
    
    
    cs_sim_data.bp_im_s{k} = bp_im
    cs_sim_data.pixelifsampled_s{k} = pixelifsampled;
    
end




cs_sim_name = strrep(raster_name, '-full.csv', 'CS-SIM-01.mat')
%%
save(fullfile(cs_sim_root, cs_sim_name), 'cs_sim_data');












