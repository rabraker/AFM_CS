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

% xdirControl = get_xdir_standard_control('const-sig');
xdirControl = get_xdir_loop_shaped_control();
%-------------- Define scan Params -------------------------
%%

npix = 512;  % image resolution.
N_prescan = 150;
width = 5;  % microns
mu_length = 0.5;  % 1000 nm. length of the horizontal mu-path.
% Unit conversions.
pix_per_micron = npix/width;
mu_pix = ceil(mu_length*pix_per_micron);

%perc_s = [0.125, 0.15, 0.25];
% freqs = [1, 2, 4, 5, 8];
freqs = [2.5];
for pp = 1:length(perc_s)
    % ************************************************
    rng(1)
    sub_sample_frac = perc_s(pp); %0.10;  % Percent of pixels to subsample.
    % ************************************************
    [pix_mask] = mu_path_mask(mu_pix, npix, npix, sub_sample_frac, false);
        
    actual_sub_sample_frac = length(find(pix_mask == 1))/npix^2;
    fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
    fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)
        
    I = ones(npix,npix)-pix_mask;
    figure(1)
    imshow(I)
            
    for k=1:length(freqs)

        raster_freq = freqs(k);
        raster_period = 1/raster_freq;

        % ************************************************
        microns_per_second = width/(raster_period/2);
        volts_per_second = microns_per_second*AFM.mic2volt_xy();
        volts_per_sample = volts_per_second * AFM.Ts;  % Ts = seconds per sample
        
        % Convert mu-path length in pixels to volts
        mu_micron = (1/pix_per_micron) * mu_pix;
        mu_volts = mu_micron * AFM.mic2volt_xy();
        mu_Nsamples = ceil(mu_volts / volts_per_sample);

        % '------------------correct by overscan for ramp ramp set ---------------
        
        N_extra = mu_overscan(xdirControl.Hyr, volts_per_sample, mu_Nsamples, N_prescan)+10;
        
        
        mpt = MuPathTraj(pix_mask, width, mu_length, microns_per_second, xdirControl.Hyr,...
            'overscan_samples', N_extra, 'pre_pad_samples', N_prescan);
        
        perc = floor(actual_sub_sample_frac*100);
        fname = mpt.get_fname();
        target_dir = sprintf('%dmicrons/parents-loop', width);
        data_root = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', target_dir);
        % data_root = fullfile(PATHS.exp(), 'z-bounce/parents')
        
        fpath_json = fullfile(data_root, fname);
        fprintf('File root:\n%s\n', data_root)
        fprintf('File name:\n%s\n', fname)
        
        mpt.write_data_json(fpath_json)
        refs = reshape(mpt.as_vector(), 3, []);
        figure(3)
        plot(refs(1, 1:1200)*512)
    end
end