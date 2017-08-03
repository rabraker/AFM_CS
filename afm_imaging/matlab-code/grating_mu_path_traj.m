

% 500 nm pitch size. The holes themselves appear to be about 1/4 of this.
% Thus, i think that if we travel for 1/2 * 500 nm = 250 nm, we should be
% able to detect wheather or not we started in a hole. The idea being that
% if no jump occures in the uz signal, we assume we are on the flat, and if
% we see a jump, we can backout which part is the hole.


width = 5;  % microns
pix = 256;

minpath = 0.25  % 250 nm.

pix_per_micron = pix/width;
mu_pix = ceil(minpath*pix_per_micron);

pixifsampled = zeros(pix, pix);

sub_sample_frac = 0.1


% ************************************************
% ************************************************

Ts = 40e-6;
Fs = 1/Ts
microns2volts = 1/5;
raster_rate = 1;  % hz
raster_period = 1/raster_rate;

microns_per_second = width/(raster_period/2)
pixels_per_second = pix_per_micron * microns_per_second
volts_per_second = microns_per_second*microns2volts;
volts_per_sample = volts_per_second * Ts  % Ts = seconds per sample

% Convert mu-path length in pixels to volts
mu_micron = (1/pix_per_micron) * mu_pix;
mu_volts = mu_micron * microns2volts;
mu_Nsamples = ceil(mu_volts / volts_per_sample)

ME = MeasEntityMu.factory([volts_per_sample, 0]);

N_mve = 1;
XR = [];
YR = [];
% ************************************************



for n=1:pix % down rows
    m = 1;
   while m < pix - mu_pix  % accros columns
       if rand(1,1) < sub_sample_frac/mu_pix
          pixifsampled(n, m:m+mu_pix) = 1;
          m = m + mu_pix;
          XR = [XR; (m /pix_per_micron) * microns2volts];
          YR = [YR; (n /pix_per_micron) * microns2volts];
       else
           m = m+1;
       end
       
       
   end
    
end

actual_sub_sample_frac = length(find(pixifsampled == 1))/pix^2;

fprintf('Desired sub sample fraction: %f\n', sub_sample_frac)
fprintf('Actual  sub sample fraction: %f\n', actual_sub_sample_frac)



% close all
I = ones(pix,pix)-pixifsampled;
I = flipud(I); % why the fuck does it start from the top???
imshow(I)
%%
meta_cell = repmat({mu_Nsamples}, 1, length(XR));
MT = MasterTrajster(XR, YR, meta_cell, MoveEntityStatic.factory(N_mve), ME);

MT.visualize_sampling;

%%


%%

mu = 1
% A X-Y points at which to take measurements. 
npts = 20;
XR = rand(npts)*1;
YR = rand(npts)*2;
% XR = zeros(1, npts);
% YR = zeros(1, npts);
% How long should we sit at eache point? 
N_me = 5000;  % Probably absurdly long. 
meta_cell = repmat({N_me}, 1, npts);

% Instantiate. 

    

