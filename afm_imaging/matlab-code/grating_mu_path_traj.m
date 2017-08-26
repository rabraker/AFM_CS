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
close all
addpath('classes')

width = 5;  % microns
pix = 256;  % image resolution.
minpath = 0.5;  % 250 nm. length of the horizontal mu-path. 
sub_sample_frac = 0.2;  % Percent of pixels to subsample. 

% Unit conversions.
pix_per_micron = pix/width;
mu_pix = ceil(minpath*pix_per_micron);

% ************************************************
% ************************************************

Ts = 40e-6;  % AFM sample rate. 
Fs = 1/Ts;
microns2volts = 1/5;
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
XR = [];
YR = [];
% ************************************************


pixifsampled = zeros(pix, pix);

for n=1:256 % down rows
    m = 1;
   while m < pix - mu_pix  % accros columns
       if rand(1,1) < sub_sample_frac/mu_pix
          pixifsampled(n, m:m+mu_pix) = 1;
          XR = [XR; ( (m - 1) / pix_per_micron) * microns2volts];
          YR = [YR; ( (n - 1) / pix_per_micron) * microns2volts];
          m = m + mu_pix;
       else
           m = m+1;
       end
   end
    
end

actual_sub_sample_frac = length(find(pixifsampled == 1))/pix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)


%
% close all
I = ones(pix,pix)-pixifsampled;
I = flipud(I); % why the fuck does it start from the top???
imshow(I)
%
meta_cell = repmat({mu_Nsamples}, 1, length(XR));
MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve), ME);

MT.visualize_sampling;

xlabel('x [v]')
ylabel('y [v]')

%%

datafile = 'C:\Users\arnold\Documents\labview\afm_imaging/data/cs-traj20-500-1Hz.csv';
if 1
    MT.write_csv(datafile)
end
   



    

