% This script reads in the full, raw data from a raster scan run. This data
% will be produced by the vi play-raster-scan.vi. That vi produces two .csv
% files. One contains the the pre-processed data that labview does in
% realtime for visualization purposed. For slow scans, using that data is
% sufficient. The vi also produces a csv file with all of the raw data. For
% faster scans, we want to use that data, so we can process it better.
% That's what this script does.
clc
clear
close all
addpath('functions')

% dat_root = '/home/arnold/gradschool/afm-cs/afm_imaging/data';
dat_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data\raster';
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data'
% datmat = csvread(fullfile(data_root, 'raster_8-1-2017_v3.csv'));

dat_name = 'raster_scan_5mic_1Hz_out_9-4-2017-01-full.csv';
dat_name = 'raster_scan_5mic_1Hz_out_9-6-2017-01-full.csv'; % Uses Dinv
dat_name = 'raster_scan_5mic_10Hz_out_9-6-2017-06-full.csv'; %uses Dinv, Ki=.02
parent_name = 'raster_scan_5mic_10Hz.csv';

dat_name = 'raster_scan_512pix_5mic_1.00e-01Hz_out_9-18-2017-03.csv'; % good
% dat_name = 'raster_scan_512pix_5mic_2.50e+00Hz_out_9-18-2017-01.csv';
% 'raster_scan_512pix_5mic_05Hz_out_9-18-2017-01.csv'
parent_name = get_parent_name(dat_name, '_out_')
meta_name = strrep(parent_name, '.csv', '.mat')

sub_dir = '5microns';
% dat_name = '';
% parent_name = 'raster_scan__256pix_20mic_1.25e-01Hz.csv';

dat_path = fullfile(dat_root, sub_dir, dat_name);
parent_path = fullfile(dat_root, sub_dir, parent_name);
meta_path = fullfile(dat_root, sub_dir,  meta_name);

tic

load(meta_path);

datmat = csvread(dat_path);
parent_dat = csvread(parent_path);


xyref = reshape(parent_dat', 2, [])';
xref = xyref(:,1);
%

% figure(1); plot(datmat(:,1));
% ax1 = gca
% figure(2); plot(datmat(:,3));
% ax3 = gca;
% figure(3); plot(datmat(:,4));
% ax4 = gca;
% linkaxes([ax1, ax3, ax4])
Ts = 40e-6;
samps_per_period = size(parent_dat,1)/2 % twice as many in here for x & y.
samps_per_line = samps_per_period/2

nperiods = 512;
pix = nperiods;
datmat = datmat([1:nperiods*samps_per_period], :);

% 
% trace_inds = get_trace_indeces(nperiods, samps_per_period);
% xdat = reshape(datmat(trace_inds,1), [], nperiods);
% ydat = reshape(datmat(trace_inds,2), [], nperiods);
% udat = reshape(datmat(trace_inds,4), [], nperiods);
% 
% udat = datmat(:,4);
% pixmat = bin_raster_fast(udat, pix, samps_per_period);
% %%
% lo = min(min(pixmat));
% hi = max(max(pixmat));
% figure(1)
% imshow(pixmat, [lo, hi])
% 
% I_fit = detrend_plane(pixmat);
% 
% figure(3)
% lo = min(min(I_fit));
% hi = max(max(I_fit));
% imshow(I_fit, [lo, hi])

% visualize tracking error.
if 0
    np = 3;
    xx = datmat(:,1);
    x_np = xx(1:np*length(xref));
    x_np = x_np - min(x_np);
    figure(200); clf; hold on
    t = [0:1:length(xref)*np-1]'*Ts;
    xref_np = repmat(xref, np, 1);
    p1 = plot(t, xref_np*volts2microns);
    p1.DisplayName = '$x_{ref}$';
    p2 = plot(t, x_np*volts2microns);
    p2.DisplayName = '$x(k)$';

    ylm = ylim;
    ylim([0, ylm(2)+.1*ylm(2)])
    ylabel('x-dir [\mu m]')
    xlabel('time [s]')
    leg1 = legend([p1, p2]);
    set(leg1, 'FontSize', 14, 'interpreter', 'latex', 'orientation', 'horizontal')
    leg1.Position=[0.6436    0.8590    0.2611    0.0640];
end
%  Now try by actually using the x-y measured data;

clc

width = 5;
% volts2micron = 50/10;
micron2pix = pix/width;
volts2pix = volts2microns * micron2pix;


pixmat2 = bin_raster_really_slow(datmat(:,[1,2,4]), pix, samps_per_period, volts2pix);

pixmat2=pixmat2(:,1:end);

pixmat2 = detrend_plane(pixmat2);
figure(5)
lo = min(min(pixmat2));
hi = max(max(pixmat2));
imshow(pixmat2, [lo, hi])

traster_recon = toc
fprintf('Total data processing time:%.3f', traster_recon)

