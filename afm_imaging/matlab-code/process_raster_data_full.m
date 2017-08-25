clc
close all
% dat_root = '/home/arnold/gradschool/afm-cs/afm_imaging/data';
dat_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data';
% data_root = 'C:\Users\arnold\Documents\labview\afm_imaging\data'
% datmat = csvread(fullfile(data_root, 'raster_8-1-2017_v3.csv'));
dat_name = 'raster-scan-8-24-2017_7Hz-full.csv';
parent_name = 'raster_traj_singleperiod_7Hz.csv';

dat_path = fullfile(dat_root, dat_name);
parent_path = fullfile(dat_root, parent_name);
datmat = csvread(dat_path);
parent_dat = csvread(parent_path);

%%

samps_per_period = size(parent_dat,1)/2; % twice as many in here for x & y.
samps_per_line = samps_per_period/2;

nperiods = 256;
pix = nperiods;
datmat = datmat([1:nperiods*samps_per_period], :);


trace_inds = get_trace_indeces(nperiods, samps_per_period);
xdat = reshape(datmat(trace_inds,1), [], nperiods);
ydat = reshape(datmat(trace_inds,2), [], nperiods);
udat = reshape(datmat(trace_inds,4), [], nperiods);

%
udat = datmat(:,4);
pixmat = bin_raster_fast(udat, pix, samps_per_period);

lo = min(min(pixmat));
hi = max(max(pixmat));
figure(1)
imshow(pixmat, [lo, hi])

I_fit = detrend_plane(pixmat);

figure(3)
lo = min(min(I_fit));
hi = max(max(I_fit));
imshow(I_fit, [lo, hi])


%% Now try by actually using the x-y measured data;



width = 5;
volts2micron = 50/10;
micron2pix = pix/width;
volts2pix = volts2micron * micron2pix;


pixmat2 = bin_raster_slow(datmat(:,[1,2,4]), pix, samps_per_period, volts2pix);



% pixmat2 = detrend_plane(pixmat2);
figure(5)
lo = min(min(pixmat2));
hi = max(max(pixmat2));
imshow(pixmat2, [lo, hi])



