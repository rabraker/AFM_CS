% The middle of our grating has a 500 nm pitch size. The holes themselves 
% appear to be about 1/4 of this. Thus, i think that if we travel for
% 1/2 * 500 nm = 250 nm, we should be able to detect wheather or not we 
% started in a hole. The idea being that if no jump occures in the u_z 
% signal, we assume we are on the flat, and if we see a jump, we can 
% backout which part is the hole. 
% 
% This will essentially come down to edge detection. If an edge is present,
% then we know the lower part of the signal was in a hole. If no edge is
% present, the entire signal was on the flat part of the grating. 
clear, clc
% close all
addpath('classes')
addpath('functions')



%-------------- Location of System Model -------------------------

plant_data = load(fullfile(PATHS.sysid(), 'x-axis_sines_infoFourierCoef_9-11-2018-01.mat'));
PLANT_init_x = ss(plant_data.modelFit.models.G_uz2stage);

%-------------- Define scan Params -------------------------
Ts = 40e-6;
TsTicks = 1600;
Ki_x = 0.01;

if 1
    width =5;  % microns
    npix =512;  % image resolution.
    mu_length = 0.5;  % 1000 nm. length of the horizontal mu-path. 
    sub_sample_frac = 0.125;  % Percent of pixels to subsample. 
end
% Unit conversions.
pix_per_micron = npix/width;
mu_pix = ceil(mu_length*pix_per_micron);

% ************************************************
% ************************************************

Ts = 40e-6;  % AFM sample rate. 
Fs = 1/Ts;
microns2volts = 1/5;  %10/50;
raster_freq = 1;  % hz
raster_period = 1/raster_freq;

microns_per_second = width/(raster_period/2);
pixels_per_second = pix_per_micron * microns_per_second;
volts_per_second = microns_per_second*microns2volts;
volts_per_sample = volts_per_second * Ts;  % Ts = seconds per sample

% Convert mu-path length in pixels to volts
mu_micron = (1/pix_per_micron) * mu_pix;
mu_volts = mu_micron * microns2volts;
mu_Nsamples = ceil(mu_volts / volts_per_sample);

ME = MeasEntityMu.factory([volts_per_sample, 0]);

N_mve = 1;

% ************************************************

[pix_mask, XR, YR] = mu_path_mask(npix, mu_pix, sub_sample_frac, width, microns2volts);

actual_sub_sample_frac = length(find(pix_mask == 1))/npix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)


I = ones(npix,npix)-pix_mask;
I = flipud(I); % why the fuck does it start from the top???
figure(1)
imshow(I)



% '------------------correct by overscan for ramp ramp set ---------------


N = mu_Nsamples;
x_rate = volts_per_sample;
G = PLANT_init_x;
D = tf([0.01, 0], [1 -1], Ts);
H = feedback(D*G, 1);
% N-1 because x0 is sample 1.
x_N =  (N-1)*x_rate;
N_extra = mu_overscan(G, D, x_rate, mu_Nsamples, 1);

mu_Nsamples_extra = N+N_extra;
meta_cell = repmat({mu_Nsamples_extra}, 1, length(XR));
%

% clc
mpt = MuPathTraj(pix_mask, width, mu_length, microns_per_second, Ts,...
  'overscan_samples', N_extra, 'pre_pad_samples', 500);
mpt.connect_mu_paths(0.064);
vec = mpt.as_vector();


perc = floor(actual_sub_sample_frac*100);
fname = sprintf('cs-traj-%dpix-%dperc-%dnm-%dmic-%.2dHz_v2.csv',npix, perc, mu_length*1000,width, raster_freq);
target_dir = sprintf('%dmicrons/parents', width);
data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', target_dir);

fpath_csv = fullfile(data_root, fname);
fprintf('File root:\n%s\n', data_root)
fprintf('File name:\n%s\n', fname)

%%
mpt.write_data(fpath_csv)
%%
% 
% MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve), ME);
% MT.visualize_sampling;
% 
% xlabel('x [v]')
% ylabel('y [v]')
% grid on

% create meta file name

target_dir = sprintf('%dmicrons/parents', width)
data_root = fullfile(getdataroot, 'cs-data', target_dir)
if exist(data_root, 'file') ~=2
    mkdir(fullfile(getdataroot, 'cs-data'), target_dir)
end

    

