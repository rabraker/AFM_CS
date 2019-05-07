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

% close all
% clear, clc
addpath('functions')
try; rmpath('functions/scanning_v0'); end
addpath('functions/scanning_v1')
addpath('functions/state_space_x')


%-------------- Location of System Model -------------------------

plants = CanonPlants.plants_ns14(9);
PLANT_init_x = plants.PLANT;

%-------------- Define scan Params -------------------------

TsTicks = 1600;
Ki_x = 0.01;

if 1
    width =5;  % microns
    npix = 512;  % image resolution.
    mu_length = 0.5;  % 1000 nm. length of the horizontal mu-path. 
    sub_sample_frac = 0.15;  % Percent of pixels to subsample. 
    raster_freq = 1;  % hz
end
N_prescan = 250;
% Unit conversions.
pix_per_micron = npix/width;
mu_pix = ceil(mu_length*pix_per_micron);

% ************************************************
% ************************************************
rng(1)

microns2volts = 1/5;  %10/50;

raster_period = 1/raster_freq;

microns_per_second = width/(raster_period/2);
volts_per_second = microns_per_second*microns2volts;
volts_per_sample = volts_per_second * AFM.Ts;  % Ts = seconds per sample

% Convert mu-path length in pixels to volts
mu_micron = (1/pix_per_micron) * mu_pix;
mu_volts = mu_micron * microns2volts;
mu_Nsamples = ceil(mu_volts / volts_per_sample);


% ************************************************

[pix_mask] = mu_path_mask(mu_pix, npix, npix, sub_sample_frac, false);

actual_sub_sample_frac = length(find(pix_mask == 1))/npix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)

Fig = mkfig(1, 4, 4);
ha = tight_subplot(1,1, 0.01, [0.01, 0.01], [0.01, 0.01]);
I = ones(npix,npix)-pix_mask;

imagesc(ha, I)
colormap('gray')

save_fig(Fig, fullfile(PATHS.thesis_root, 'figures/implementation_0/mu_path_mask_512'), false)

